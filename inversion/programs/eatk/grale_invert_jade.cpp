#include "log.h"
#include "inputoutput.h"
#include "inversioncommunicatornewga.h"
#include "lensgaindividual.h"
#include "lensgastopcriterion.h"
#include "lensfitnessgeneral.h"
#include "lensgaconvergenceparameters.h"
#include "lensgafitnesscomparison.h"
#include "lensdemutationcrossover.h"
#include "lensdeevolver.h"
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

	auto mut = make_shared<grale::LensDEMutation>();
	auto cross = make_shared<grale::LensDECrossover>(rng, genomeCalculator->allowNegativeValues()); // TODO: same

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

	grale::LensDEEvolver jade(rng, mut, cross, comparison); // TODO: make other JADE parameters configurable?

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
