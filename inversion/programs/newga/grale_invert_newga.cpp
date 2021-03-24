#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include "lensinversiongafactorycommon.h"
#include "inversioncommunicatornewga.h"
#include "randomnumbergenerator.h"
#include "constants.h"
#include "lensgaindividual.h"
#include "lensgagenomemutation.h"
#include "lensgagenomecrossover.h"
#include "lensgafitnesscomparison.h"
#include "lensgasingleobjectivecrossover.h"
#include "lensgamultiobjectivecrossover.h"
#include "galensmodule.h"
#include "lensgastopcriterion.h"
#include <serut/memoryserializer.h>
#include <mogal/geneticalgorithm.h>
#include <mogal2/geneticalgorithm.h>
#include <mogal2/randomnumbergenerator.h>
#include <mogal2/singlethreadedpopulationfitnesscalculation.h>
#include <mogal2/multithreadedpopulationfitnesscalculation.h>
#include <mogal2/stopcriterion.h>
#include <mogal2/simplesortedpopulation.h>
#include <serut/vectorserializer.h>
#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#endif // !WIN32

#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <limits>
#include <chrono>

using namespace std;
using namespace serut;
using namespace errut;

// Wrap the GSL based RNG for now
class RNG : public mogal2::RandomNumberGenerator
{
public:
	RNG() { }
	~RNG() { }
    double getRandomDouble() override { return m_rng.pickRandomNumber(); }
    float getRandomFloat() override { return (float)m_rng.pickRandomNumber(); }
	uint32_t getRandomUint32() override
	{
		cerr << "getRandomUint32 NOT IMPLEMENTED" << endl;
		exit(-1);
		return 0;
	}
private:
	grale::RandomNumberGenerator m_rng;
};

class LensGAMultiObjectiveCrossover : public mogal2::PopulationCrossover
{

};

class LensFitnessCalculation : public mogal2::GenomeFitnessCalculation
{
public:
	LensFitnessCalculation(grale::LensInversionGAFactoryCommon &factory) : m_pFactory(&factory) { }
	LensFitnessCalculation(const shared_ptr<grale::LensInversionGAFactoryCommon> &factory)
	{
		m_spFactory = factory; // keep a reference to the lib
		m_pFactory = factory.get();
	}

	~LensFitnessCalculation() { }

	errut::bool_t calculate(const mogal2::Genome &genome, mogal2::Fitness &fitness)
	{
		const grale::LensGAGenome &g = dynamic_cast<const grale::LensGAGenome &>(genome);
		grale::LensGAFitness &f = dynamic_cast<grale::LensGAFitness &>(fitness);
		
		// Anything that needs to get communicated back needs to be in the fitness, so
		// we'll store the scalefactor in the fitness and transfer it to the genome in a
		// later stage
		if (!m_pFactory->calculateFitness(g.m_weights, g.m_sheets, f.m_scaleFactor, f.m_fitnesses.data()))
			return "Error calculating fitness: " + m_pFactory->getErrorString();


		// TODO: for testing
		// usleep(10000);

		return true;
	}
private:
	grale::LensInversionGAFactoryCommon *m_pFactory;
	shared_ptr<grale::LensInversionGAFactoryCommon> m_spFactory;
};

class Timer {
public:
    Timer() {
        start();
    }

    void start() {
        beg = std::chrono::steady_clock::now();
    }

    void stop() {
        end = std::chrono::steady_clock::now();
    }

    double duration() {
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - beg).count();
    }
private:
    std::chrono::time_point<std::chrono::steady_clock> beg;
    std::chrono::time_point<std::chrono::steady_clock> end;
};


class MyGA : public mogal2::GeneticAlgorithm
{
public:
	MyGA() { }
	~MyGA()
	{
		double dtSum = 0;
		double dt2Sum = 0;
		for (auto dt : m_intervals)
		{
			dtSum += dt;
			dt2Sum += dt*dt;
		}

		double avg = dtSum/m_intervals.size();
		double stddev = SQRT(dt2Sum/m_intervals.size() - avg*avg);
		cerr << "Avg = " << avg/1e6 << " ms, stddev = " << stddev/1e6 << " ms" << endl;
	}

	Timer m_timer;
	vector <double> m_intervals;

	errut::bool_t onBeforeFitnessCalculation(size_t generation, std::shared_ptr<mogal2::Population> &population)
	{
		m_timer.start();
	// 	cout << "# Generation " << generation << ", before calculation: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, std::shared_ptr<mogal2::Population> &population)
	{
		m_timer.stop();
		m_intervals.push_back(m_timer.duration());
	// 	cout << "# Generation " << generation << ", calculated: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
		return true;
	}
};

class PreferredIndividualSelector
{
public:
	PreferredIndividualSelector() { }
	virtual ~PreferredIndividualSelector() { }
	virtual errut::bool_t select(const std::vector<std::shared_ptr<mogal2::Individual>> &best,
								 std::shared_ptr<mogal2::Individual> &selected)
	{
		return "Not implemented";
	}
};

class SubsequentBestIndividualSelector : public PreferredIndividualSelector
{
public:
	SubsequentBestIndividualSelector(size_t numObjectives,
									 const std::shared_ptr<mogal2::FitnessComparison> &fitComp)
		: m_numObjectives(numObjectives), m_cmp(fitComp) { }
	~SubsequentBestIndividualSelector() { }
	
	errut::bool_t select(const std::vector<std::shared_ptr<mogal2::Individual>> &bestIndividuals,
								 std::shared_ptr<mogal2::Individual> &selected) override
	{
		vector<shared_ptr<mogal2::Individual>> genomes = bestIndividuals;
		vector<shared_ptr<mogal2::Individual>> genomes2;
		shared_ptr<mogal2::Individual> best;

		if (genomes.size() == 0)
			return "No best individuals to select one from";

		for (size_t comp = 0 ; comp < m_numObjectives ; comp++)
		{
			best = genomes[0];

			for (auto &i : genomes)
			{
				if (m_cmp->isFitterThan(i->fitnessRef(), best->fitnessRef(), comp))
					best = i;
			}

			genomes2.clear();

			// Ok, now we know a genome that has the lowest 'comp' fitness, let's see it there
			// are others which perfom equally well
			for (auto &i : genomes)
			{
				// if 'best' is not fitter than 'i' it must have the same fitness with respect to this component
				if (!m_cmp->isFitterThan(best->fitnessRef(), i->fitnessRef(), comp))
					genomes2.push_back(i);
			}

			swap(genomes, genomes2);
			genomes2.clear();
		}
		selected = best;
		return true;
	}
private:
	size_t m_numObjectives;
	std::shared_ptr<mogal2::FitnessComparison> m_cmp;
};

// TODO: rename this, is from copy-paste
class SingleCoreCommunicator : public InversionCommunicator
{
public:
	SingleCoreCommunicator() { }
	~SingleCoreCommunicator() { }
protected:
	string getVersionInfo() const override { return "Single core genetic algorithm engine"; }
	
	bool hasGeneticAlgorithm() const override { return true; }
	void getAllBestGenomes(std::vector<std::shared_ptr<mogal2::Individual>> &bestGenomes) override
	{
		bestGenomes = m_best;
	}

	shared_ptr<mogal2::Individual> getPreferredBestGenome() override
	{
		if (m_best.size() == 0)
			return nullptr;
		if (!m_selector.get())
			return nullptr;
		
		bool_t r;
		shared_ptr<mogal2::Individual> selected;
		if (!(r = m_selector->select(m_best, selected)))
		{
			cerr << "Error in preferred individual selection: " << r.getErrorString() << endl;
			return nullptr;
		}

		return selected;
	}

	bool_t runGA(int popSize, mogal::GAFactory &factory, mogal::GeneticAlgorithmParams &params,
	             const std::string &moduleDir, const std::string &moduleFile,
	             const std::vector<uint8_t> &factoryParamBytes) override
	{
		bool_t r;
		grale::LensInversionGAFactoryCommon &gaFactory = dynamic_cast<grale::LensInversionGAFactoryCommon&>(factory);
		size_t numObjectives = gaFactory.getNumberOfFitnessComponents();

		m_selector = make_shared<SubsequentBestIndividualSelector>(numObjectives,
								make_shared<grale::LensGAFitnessComparison>());

		shared_ptr<RNG> rng = make_shared<RNG>();
		MyGA ga;

		auto mutation = make_shared<grale::LensGAGenomeMutation>(rng, 
					   gaFactory.getChanceMultiplier(),
					   gaFactory.allowNegativeValues(),
					   gaFactory.getMutationAmplitude(),
					   gaFactory.useAbsoluteMutation());

		grale::LensGAIndividualCreation creation(rng, gaFactory.getNumberOfBasisFunctions(), gaFactory.getNumberOfSheets(),
		                  gaFactory.allowNegativeValues(), gaFactory.getNumberOfFitnessComponents());

		shared_ptr<grale::LensGACrossoverBase> cross;
		if (numObjectives == 1)
			cross = make_shared<grale::LensGASingleObjectiveCrossover>(params.getBeta(),
					      params.useElitism(),
						  params.alwaysIncludeBestGenome(),
						  params.getCrossOverRate(),
						  rng,
						  gaFactory.allowNegativeValues(),
						  mutation);
		else
			cross = make_shared<grale::LensGAMultiObjectiveCrossover>(params.getBeta(),
					      params.useElitism(),
						  params.alwaysIncludeBestGenome(),
						  params.getCrossOverRate(),
						  rng,
						  gaFactory.allowNegativeValues(),
						  mutation,
						  numObjectives);

		size_t numThreads = 1;
		if (getenv("NUMTHREADS"))
			numThreads = stoi(getenv("NUMTHREADS"));
		cerr << "Using " << numThreads << " threads " << endl;

		vector<shared_ptr<mogal2::GenomeFitnessCalculation>> genomeFitnessCalculators = { make_shared<LensFitnessCalculation>(gaFactory) };
		grale::GALensModule module;

		if (!module.open(moduleDir, moduleFile))
			return "Couldn't open module: " + module.getErrorString();

		for (size_t i = 2 ; i < numThreads ; i++)
		{
			auto newFac = shared_ptr<grale::LensInversionGAFactoryCommon>(dynamic_cast<grale::LensInversionGAFactoryCommon*>(module.createFactoryInstance()));
			auto facParams = shared_ptr<mogal::GAFactoryParams>(newFac->createParamsInstance());

			serut::VectorSerializer ser(factoryParamBytes);
			if (!facParams->read(ser))
				return "Coudln't deserialize factory parameters: " + facParams->getErrorString();

			if (!newFac->init(facParams.get()))
				return "Unable to init factory copy: " + newFac->getErrorString();
			
			genomeFitnessCalculators.push_back(make_shared<LensFitnessCalculation>(newFac));
		}
		
		shared_ptr<mogal2::PopulationFitnessCalculation> calc;
		if (numThreads <= 1)
			calc = make_shared<mogal2::SingleThreadedPopulationFitnessCalculation>(genomeFitnessCalculators[0]);
		else
		{
			auto mpCalc = make_shared<mogal2::MultiThreadedPopulationFitnessCalculation>();
			if (!(r = mpCalc->initThreadPool(genomeFitnessCalculators)))
				return "Unable to initialize threads: " + r.getErrorString();
			calc = mpCalc;
		}
		
		grale::LensGAStropCriterion stop(gaFactory.getMaximumNumberOfGenerations(), mutation);
		if (!(r = stop.initialize(gaFactory.getFitnessObject())))
			return "Error initializing convergence checker: " + r.getErrorString();

		if (!(r = ga.run(creation, *cross, *calc, stop, popSize)))
			return "Error running GA: " + r.getErrorString();

		m_best = cross->getBestIndividuals();
		cout << "Best: " << endl;
		for (auto &b: m_best)
			cout << b->fitness()->toString() << endl;

		return true;
	}
	
	std::vector<std::shared_ptr<mogal2::Individual>> m_best;
	shared_ptr<SubsequentBestIndividualSelector> m_selector;
};

int main(int argc, char *argv[])
{
	grale::LOG.init(argv[0]);
#ifndef WIN32
	// This was added to be able to launch instance in gdb, in which case
	// stdin/stdout can no longer be used for communication with the
	// executable (needed to control gdb).
	if (argc == 3)
	{
		string inputPipeName(argv[1]);
		string outputPipeName(argv[2]);
		cerr << "Opening " << outputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdOutFileDescriptor = open(outputPipeName.c_str(), O_WRONLY); 
		cerr << "Opening " << inputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = open(inputPipeName.c_str(), O_RDWR); 
	}
#endif // WIN32

	SingleCoreCommunicator comm;

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
