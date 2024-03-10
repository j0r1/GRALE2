#include "log.h"
#include "inputoutput.h"
#include "inversioncommunicatornewga.h"
#include "lensgaindividual.h"
#include "lensgastopcriterion.h"
#include "lensfitnessgeneral.h"
#include "lensgaconvergenceparameters.h"
#include "lensgafitnesscomparison.h"
#include "randomnumbergenerator.h"

#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#endif // !WIN32
#include <memory>
#include <iostream>
#include <limits>

#include <serut/memoryserializer.h>
#include <serut/vectorserializer.h>
#include <eatk/evolutionaryalgorithm.h>
#include <eatk/jadeevolver.h>
#include <eatk/singlethreadedpopulationfitnesscalculation.h>
#include <eatk/multithreadedpopulationfitnesscalculation.h>

using namespace std;
using namespace serut;
using namespace errut;

// TODO: copy from newgacommunicatorbase.h, move to separate file
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

// TODO: copy from newgacommunicatorbase.h, move to separate file
class MyGA : public eatk::EvolutionaryAlgorithm
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
		std::cerr << "Avg = " << avg/1e6 << " ms, stddev = " << stddev/1e6 << " ms" << std::endl;
	}

	Timer m_timer;
	std::vector <double> m_intervals;

	errut::bool_t onBeforeFitnessCalculation(size_t generation, const std::shared_ptr<eatk::Population> &population) override
	{
		m_timer.start();
	// 	cout << "# Generation " << generation << ", before calculation: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, const std::shared_ptr<eatk::Population> &population) override
	{
		m_timer.stop();
		m_intervals.push_back(m_timer.duration());
	 	//cerr << "# Generation " << generation << ", calculated: " << endl;
	 	//for (auto &i : population->individuals())
	 	//	cerr << i->fitness()->toString() << " "; // endl;
		//cerr << endl;
		return true;
	}

	errut::bool_t onBeforeFitnessCalculation(size_t generation, const std::vector<std::shared_ptr<eatk::Population>> &populations) override
	{
		m_timer.start();
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, const std::vector<std::shared_ptr<eatk::Population>> &populations) override
	{
		m_timer.stop();
		m_intervals.push_back(m_timer.duration());
		return true;
	}
};

class DEMutation : public eatk::DifferentialEvolutionMutation
{
public:
	DEMutation() { }
	~DEMutation() { }

	errut::bool_t check(const eatk::Genome &g) override
	{
		if (!dynamic_cast<const grale::LensGAGenome*>(&g))
			return "Genome is of wrong type";
		return true;
	}

	std::shared_ptr<eatk::Genome> mutate(const std::vector<const eatk::Genome*> &genomes, const std::vector<double> &weights) override
	{
		assert(genomes.size() > 0);
		assert(genomes.size() == weights.size());

		std::shared_ptr<eatk::Genome> result = genomes[0]->createCopy();
		grale::LensGAGenome &g0 = static_cast<grale::LensGAGenome&>(*result);

		// Initialize everything to zero, so we can safely add things
		g0.m_weights.assign(g0.m_weights.size(), 0);
		g0.m_sheets.assign(g0.m_sheets.size(), 0);
	
		// We can't set this to NaN, since it will be used in the next crossover step
		// But we've incorporated the scalefactors, so this should just be 1
		g0.m_scaleFactor = 1.0;

		for (size_t j = 0 ; j < genomes.size() ; j++)
		{
			float f = (float)weights[j];
			const grale::LensGAGenome &g1 = static_cast<const grale::LensGAGenome&>(*genomes[j]);

			// Add the weights, need to take the scale factor into account as well!
			assert(g0.m_weights.size() == g1.m_weights.size());
			for (size_t i = 0 ; i < g0.m_weights.size() ; i++)
				g0.m_weights[i] += f * g1.m_weights[i] * g1.m_scaleFactor;

			// Add the sheet values, no scale factor here
			assert(g0.m_sheets.size() == g1.m_sheets.size());
			for (size_t i = 0 ; i < g0.m_sheets.size() ; i++)
				g0.m_sheets[i] += f * g1.m_sheets[i];
		}

		return result;
	}
};

class DECrossover : public eatk::DifferentialEvolutionCrossover
{
public:
	DECrossover(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, bool allowNeg) : m_rng(rng), m_allowNegative(allowNeg) { }
	~DECrossover() { }

	errut::bool_t check(const eatk::Genome &g) override
	{
		if (!dynamic_cast<const grale::LensGAGenome*>(&g))
			return "Genome is of wrong type";
		return true;
	}

	errut::bool_t crossover(double CR, eatk::Genome &mutantDest0, const eatk::Genome &origVector0) override
	{
		grale::LensGAGenome &g0 = static_cast<grale::LensGAGenome&>(mutantDest0);
		const grale::LensGAGenome &g1 = static_cast<const grale::LensGAGenome&>(origVector0);

		size_t numWeigths = g0.m_weights.size();
		size_t numSheets = g0.m_sheets.size();
		size_t totalSize = numWeigths + numSheets;
		size_t rndIdx = ((size_t)m_rng->getRandomUint32())%(totalSize);

		auto getValue = [numWeigths,numSheets](const grale::LensGAGenome &g, size_t i)
		{
			if (i < numWeigths)
				return g.m_weights[i] * g.m_scaleFactor; // Take scale factor into account!
			
			i -= numWeigths;
			assert(i < numSheets);
			return g.m_sheets[i]; // no scale factor here
		};

		// Note that no scale factor will be used here
		auto setValue = [numWeigths,numSheets](grale::LensGAGenome &g, size_t i, float value)
		{
			if (i < numWeigths)
				g.m_weights[i] = value;
			else
			{
				i -= numWeigths;
				assert(i < numSheets);
				g.m_sheets[i] = value;
			}
		};

		for (size_t i = 0 ; i < totalSize ; i++)
		{
			float val = std::numeric_limits<float>::quiet_NaN();

			// We make sure to get the value in either case, so we get something that incorporated
			// the scale factor; we can then store the value again
			if (i != rndIdx && m_rng->getRandomDouble() > CR)
				val = getValue(g1, i);
			else
				val = getValue(g0, i);

			setValue(g0, i, val);
		}

		// Reset the scale factor
		g0.m_scaleFactor = std::numeric_limits<float>::quiet_NaN();

		if (!m_allowNegative) // Enforce bounds on the weights
		{
			// g1 should already be within bounds, we won't check this (only in assert)!
			for (size_t i = 0 ; i < numWeigths ; i++)
			{
				if (g0.m_weights[i] < 0)
					g0.m_weights[i] = g1.m_weights[i]/2.0f;

				assert(g0.m_weights[i] >= 0);
				assert(g1.m_weights[i] >= 0);
			}
		}

		// Always enforce bounds on the sheets - TODO: what is done in the normal GA?
		for (size_t i = 0 ; i < numSheets ; i++)
		{
			// g1 should already be within bounds, we won't check this (only in assert)!
			if (g0.m_sheets[i] < 0)
				g0.m_sheets[i] = g1.m_sheets[i]/2.0f;

			assert(g0.m_sheets[i] >= 0);
			assert(g1.m_sheets[i] >= 0);
		}

		return true;
	}

private:
	std::shared_ptr<eatk::RandomNumberGenerator> m_rng;
	const bool m_allowNegative;
};

class JADEInversionCommunicator : public InversionCommunicator
{
public:
	JADEInversionCommunicator(int numThreads) : m_numThreads(numThreads)
	{
		cerr << "Using " << m_numThreads << " threads" << endl;
	}
	~JADEInversionCommunicator() { }

	string getVersionInfo() const { return "EATk based JADE algorithm"; }

	virtual void calculatorCleanup() { } // TODO: for the future, for MPI

	bool_t runGA(int popSize, const std::string &lensFitnessObjectType,
						 const std::string &calculatorType,
	                     grale::LensGACalculatorFactory &calcFactory, 
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
						 const std::vector<uint8_t> &factoryParamBytes,
						 const grale::GAParameters &params,
						 const grale::LensGAConvergenceParameters &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams);
protected:
	// TODO: copied this
	bool_t getCalculator(const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									grale::LensGAIndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) /* override */
	{
		return getMultiThreadedPopulationCalculator(m_numThreads, lensFitnessObjectType, calculatorType,
				                                    calcFactory, genomeCalculator, factoryParamBytes, calc);
	}

	bool hasGeneticAlgorithm() const override { return true; }
	virtual void getAllBestGenomes(std::vector<std::shared_ptr<eatk::Individual>> &bestGenomes) { bestGenomes = m_best; }
	virtual std::shared_ptr<eatk::Individual> getPreferredBestGenome()
	{ 
		if (m_best.size() == 0)
			return nullptr;
		return m_best[0]; // TODO: for now, only single objective, so a single best one
	}

private:
	// TODO: just copied from newgacommunicatorbase.h, use common code for this
	static bool_t getMultiThreadedPopulationCalculator(size_t numThreads, const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc)
	{
		if (numThreads <= 1)
		{
			calc = std::make_shared<eatk::SingleThreadedPopulationFitnessCalculation>(genomeCalculator);
			return true;
		}

		std::vector<std::shared_ptr<eatk::GenomeFitnessCalculation>> genomeFitnessCalculators;
		genomeFitnessCalculators.resize(numThreads);
		genomeFitnessCalculators[0] = genomeCalculator;

		auto createAndInitCalculator = [&genomeFitnessCalculators,&lensFitnessObjectType, &calcFactory, &factoryParamBytes](size_t idx) -> bool_t
		{
			std::unique_ptr<grale::LensFitnessObject> fitObj = grale::LensFitnessObjectRegistry::instance().createFitnessObject(lensFitnessObjectType);
			if (!fitObj.get())
				return "No fitness object with name '" + lensFitnessObjectType + "' is known";

			bool_t r;
			// TODO: should we do this beforehand to make really really sure everything is thread safe?
			auto calculatorParams = calcFactory.createParametersInstance();
			if (!(r = InversionCommunicator::loadFromBytes(*calculatorParams, factoryParamBytes)))
				return "Can't load calculator parameters from received data: " + r.getErrorString();
		
			std::shared_ptr<grale::LensGAGenomeCalculator> calculatorInstance = calcFactory.createCalculatorInstance(move(fitObj));
			if (!(r = calculatorInstance->init(*calculatorParams)))
				return "Unable to initialize calculator: " + r.getErrorString();

			assert(idx < genomeFitnessCalculators.size());
			assert(!genomeFitnessCalculators[idx].get());
			genomeFitnessCalculators[idx] = calculatorInstance;
			return true;
		};

		std::vector<std::thread> calcInitThreads(numThreads-1); // The first one is already initialized, do the rest in parallel
		std::vector<bool_t> errors(calcInitThreads.size());

		auto createAndInitCalculatorWrapper = [&createAndInitCalculator,&errors](size_t idxMinusOne)
		{
			errors[idxMinusOne] = createAndInitCalculator(idxMinusOne+1);
		};

		for (size_t i = 0 ; i < calcInitThreads.size() ; i++)
			calcInitThreads[i] = std::thread(createAndInitCalculatorWrapper, i);

		for (auto &t : calcInitThreads)
			t.join();

		for (auto &r : errors)
			if (!r)
				return r;

		bool_t r;
		auto mpCalc = std::make_shared<eatk::MultiThreadedPopulationFitnessCalculation>();
		if (!(r = mpCalc->initThreadPool(genomeFitnessCalculators)))
			return "Unable to initialize threads: " + r.getErrorString();
		calc = mpCalc;

		return true;
	}

	int m_numThreads;
	std::vector<std::shared_ptr<eatk::Individual>> m_best;
};

// TODO: copy from newgacommunicatorbase.h, put in shared file
class Stop : public grale::LensGAStopCriterion
{
public:
	Stop(const std::shared_ptr<grale::LensGAGenomeMutation> &mutation, int popId = -1)
		: grale::LensGAStopCriterion(mutation), m_popId(popId) { }
protected:
	void onReport(const std::string &s)	const override
	{
		if (m_popId < 0)
			WriteLineStdout("GAMESSAGESTR:" + s);
		else
			WriteLineStdout("GAMESSAGESTR: P(" + std::to_string(m_popId) + "):" + s);
	}
private:
	int m_popId;
};

class LensGAJADEEvolver : public eatk::JADEEvolver
{
public:
	LensGAJADEEvolver(
		const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
		const std::shared_ptr<eatk::DifferentialEvolutionMutation> &mut,
		const std::shared_ptr<eatk::DifferentialEvolutionCrossover> &cross,
		const std::shared_ptr<eatk::FitnessComparison> &fitComp) : JADEEvolver(rng, mut, cross, fitComp) { }

	// We need to override this function to copy the calculated scale factors to the genomes
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
	{
		for (auto &i : population->individuals())
		{
			auto &g = static_cast<grale::LensGAGenome &>(i->genomeRef());
			auto &f = static_cast<grale::LensGAFitness &>(i->fitnessRef());
			g.m_scaleFactor = f.m_scaleFactor;
		}
		return eatk::JADEEvolver::createNewPopulation(generation, population, targetPopulationSize);
	}
};

bool_t JADEInversionCommunicator::runGA(int popSize, const std::string &lensFitnessObjectType,
						 const std::string &calculatorType,
	                     grale::LensGACalculatorFactory &calcFactory, 
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
						 const std::vector<uint8_t> &factoryParamBytes,
						 const grale::GAParameters &params,
						 const grale::LensGAConvergenceParameters &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams)
{
	errut::bool_t r;

	int numObj = genomeCalculator->getNumberOfObjectives();
	if (numObj != 1)
		return "JADE currently only works with one objective, but there are " + to_string(numObj);

	std::shared_ptr<grale::RandomNumberGenerator> rng = std::make_shared<grale::RandomNumberGenerator>();

	// TODO: At the moment we need this for the stop criterion
	std::shared_ptr<grale::LensGAGenomeMutation> mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 
					   1.0, // chance multiplier; has always been set to one
					   genomeCalculator->allowNegativeValues(),
					   0, // mutation amplitude, will be in the stop criterion
					   true); // absolute or small mutation, will be set in the stop criterion

	// TODO: this stop criterion where the mutation amplitude is changed doesn't really
	//       make sense (not the kind of mutation in DE)

	Stop stop(mutation);
	if (!(r = stop.initialize(numObj, convParams)))
		return "Can't initialize stop criterion: " + r.getErrorString();

	auto mut = make_shared<DEMutation>();
	auto cross = make_shared<DECrossover>(rng, genomeCalculator->allowNegativeValues()); // TODO: same

	auto comparison = std::make_shared<grale::LensGAFitnessComparison>();

	grale::LensGAIndividualCreation creation(rng, 
					  genomeCalculator->getNumberOfBasisFunctions(),
					  genomeCalculator->getNumberOfSheets(),
					  genomeCalculator->allowNegativeValues(),
					  genomeCalculator->getNumberOfObjectives());

	// After this is created, calculatorCleanup() should be called as well (for MPI at the moment)
	std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
	if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
							factoryParamBytes, creation, calc)))
		return "Can't get calculator: " + r.getErrorString();

	LensGAJADEEvolver jade(rng, mut, cross, comparison); // TODO: make other JADE parameters configurable?

	MyGA ga;
	if (!(r = ga.run(creation, jade, *calc, stop, popSize, popSize, popSize*2)))
	{
		calculatorCleanup();
		return "Error running GA: " + r.getErrorString();
	}

	m_best = jade.getBestIndividuals();

	return true;
}

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

	size_t numThreads = 1;
	if (getenv("GRALE_NUMTHREADS"))
		numThreads = stoi(getenv("GRALE_NUMTHREADS"));

	JADEInversionCommunicator comm(numThreads);

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
