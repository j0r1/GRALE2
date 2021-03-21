#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include "lensinversiongafactorycommon.h"
#include "inversioncommunicator.h"
#include "randomnumbergenerator.h"
#include "constants.h"
#include "lensgaindividual.h"
#include "lensgagenomemutation.h"
#include "lensgagenomecrossover.h"
#include "lensgafitnesscomparison.h"
#include "lensgasingleobjectivecrossover.h"
#include "lensfitnessgeneral.h"
#include "multifitnesshistory.h"
#include "galensmodule.h"
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
#endif // !WIN32

#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <limits>

using namespace std;
using namespace serut;
using namespace errut;

// Wrap the GSL based RNG for now
class RNG : public mogal2::RandomNumberGenerator
{
public:
	RNG(const mogal::RandomNumberGenerator *pRng) : m_pRng(pRng) { }
	~RNG() { }
    double getRandomDouble() override { return m_pRng->pickRandomNumber(); }
    float getRandomFloat() override { return (float)m_pRng->pickRandomNumber(); }
	uint32_t getRandomUint32() override
	{
		cerr << "getRandomUint32 NOT IMPLEMENTED" << endl;
		exit(-1);
		return 0;
	}
private:
	const mogal::RandomNumberGenerator *m_pRng;
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
		return true;
	}
private:
	grale::LensInversionGAFactoryCommon *m_pFactory;
	shared_ptr<grale::LensInversionGAFactoryCommon> m_spFactory;
};

class LensGAStropCriterion : public mogal2::StopCriterion
{
public:
	LensGAStropCriterion(size_t maxGenerations,
						 const std::shared_ptr<grale::LensGAGenomeMutation> &mutation)
		: m_maxGenerations(maxGenerations), m_mutation(mutation)
	{ 
		m_numObjectives = 0; // means not initialized
	}

	errut::bool_t initialize(const grale::LensFitnessObject &fitnessObject)
	{
		if (m_numObjectives > 0)
			return "Already initialized";

		auto pFitnessObject = dynamic_cast<const grale::LensFitnessGeneral*>(&fitnessObject);
		if (!pFitnessObject)
			return "Not a LensFitnessGeneral object";

		m_numObjectives = fitnessObject.getNumberOfFitnessComponents();

		int histSize = pFitnessObject->getConvergenceHistorySize();
		m_fitnessConvergenceFactors = pFitnessObject->getConvergenceFactors();
		m_mutationSizes = pFitnessObject->getConvergenceSmallMutationSizes();

		if (m_fitnessConvergenceFactors.size() != m_mutationSizes.size() || m_fitnessConvergenceFactors.size() < 1)
			return "Unexpected: invalid convergence or mutation settings (should have been checked before)";
		
		m_pFitnessHistory = make_unique<grale::MultiFitnessHistory>(m_numObjectives, histSize, 0); // Convergence factor will be reset in the next subroutine
		m_convergenceFactorPos = 0;
		updateMutationAndConvergenceInfo();

		return true;
	}

	void updateMutationAndConvergenceInfo()
	{
		float mutationSize = m_mutationSizes[m_convergenceFactorPos];
		bool useAbsoluteMutation = (mutationSize < 0)?true:false;

		m_mutation->setAbsoluteMutation(useAbsoluteMutation);
		m_mutation->setMutationAmplitude(mutationSize);
		
		m_pFitnessHistory->reset(m_fitnessConvergenceFactors[m_convergenceFactorPos]);
		m_convergenceFactorPos++;

		if (!useAbsoluteMutation)
			cerr << "DEBUG: Setting small mutation with size " << to_string(mutationSize) << endl;
		else
			cerr << "DEBUG: Setting large mutations" << endl;
	}

	errut::bool_t analyze(const std::vector<std::shared_ptr<mogal2::Individual>> &currentBest, size_t generationNumber, bool &shouldStop) override
	{
		if (m_numObjectives == 0)
			return "Not initialized";

		assert(m_pFitnessHistory.get());

		for (int i = 0 ; i < m_numObjectives ; i++)
		{
			for (int j = 0 ; j < currentBest.size() ; j++)
			{
				const grale::LensGAFitness &f = static_cast<const grale::LensGAFitness&>(currentBest[i]->fitnessRef());

				assert(i < (int)f.m_fitnesses.size());
				m_pFitnessHistory->processValue(i, f.m_fitnesses[i]);
			}
		}

		std::stringstream ss;
		// Logging generation-1 for comparison with old code
		ss << "Generation " << (generationNumber-1) << ": " << m_pFitnessHistory->getDebugInfo();
		cerr << ss.str() << endl;

		if (m_pFitnessHistory->isFinished())
		{
			if (m_convergenceFactorPos == m_fitnessConvergenceFactors.size())
				shouldStop = true;
			else
				updateMutationAndConvergenceInfo();
		}

		m_pFitnessHistory->advance();

		if (generationNumber > m_maxGenerations)
			shouldStop = true;

		return true;
	}
private:
	size_t m_maxGenerations;
	std::shared_ptr<grale::LensGAGenomeMutation> m_mutation;
	size_t m_numObjectives;
	std::vector<double> m_fitnessConvergenceFactors;
	std::vector<double> m_mutationSizes;
	std::unique_ptr<grale::MultiFitnessHistory> m_pFitnessHistory;
	int m_convergenceFactorPos;
};

class MyGA : public mogal2::GeneticAlgorithm
{
public:
	MyGA() { }
	~MyGA() { }

	// errut::bool_t onBeforeFitnessCalculation(size_t generation, std::shared_ptr<mogal2::Population> &population)
	// {
	// 	cout << "# Generation " << generation << ", before calculation: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
	// 	return true;
	// }

    // errut::bool_t onFitnessCalculated(size_t generation, std::shared_ptr<mogal2::Population> &population)
	// {
	// 	cout << "# Generation " << generation << ", calculated: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
	// 	return true;
	// }
};

// TODO: rename this, is from copy-paste
class SingleCoreCommunicator : public InversionCommunicator
{
public:
	SingleCoreCommunicator() { }
	~SingleCoreCommunicator() { }
protected:
	string getVersionInfo() const override { return "Single core genetic algorithm engine"; }

	bool_t runGA(int popSize, mogal::GAFactory &factory, mogal::GeneticAlgorithmParams &params,
	             const std::string &moduleDir, const std::string &moduleFile,
	             const std::vector<uint8_t> &factoryParamBytes) override
	{
		bool_t r;
		grale::LensInversionGAFactoryCommon &gaFactory = dynamic_cast<grale::LensInversionGAFactoryCommon&>(factory);
		if (gaFactory.getNumberOfFitnessComponents() != 1)
			return "Currently only for one fitness component";

		shared_ptr<RNG> rng = make_shared<RNG>(gaFactory.getRandomNumberGenerator());
		MyGA ga;

		auto mutation = make_shared<grale::LensGAGenomeMutation>(rng, 
					   gaFactory.getChanceMultiplier(),
					   gaFactory.allowNegativeValues(),
					   gaFactory.getMutationAmplitude(),
					   gaFactory.useAbsoluteMutation());

		grale::LensGAIndividualCreation creation(rng, gaFactory.getNumberOfBasisFunctions(), gaFactory.getNumberOfSheets(),
		                  gaFactory.allowNegativeValues(), gaFactory.getNumberOfFitnessComponents());

		grale::LensGASingleFitnessCrossover cross(params.getBeta(),
					      params.useElitism(),
						  params.alwaysIncludeBestGenome(),
						  params.getCrossOverRate(),
						  rng,
						  gaFactory.allowNegativeValues(),
						  mutation);

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
		
		LensGAStropCriterion stop(gaFactory.getMaximumNumberOfGenerations(), mutation);
		if (!(r = stop.initialize(gaFactory.getFitnessObject())))
			return "Error initializing convergence checker: " + r.getErrorString();

		if (!(r = ga.run(creation, cross, *calc, stop, popSize)))
			return "Error running GA: " + r.getErrorString();

		auto &bestSolutions = cross.getBestIndividuals();
		cout << "Best: " << endl;
		for (auto &b: bestSolutions)
			cout << b->fitness()->toString() << endl;

		return true;
	}
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
