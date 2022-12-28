#pragma once

#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactorycommon.h"
#include "inversioncommunicatornewga.h"
#include "constants.h"
#include "lensgaindividual.h"
#include "lensgagenomemutation.h"
#include "lensgagenomecrossover.h"
#include "lensgafitnesscomparison.h"
#include "lensgasingleobjectivecrossover.h"
#include "lensgamultiobjectivecrossover.h"
#include "lensgastopcriterion.h"
#include "lensfitnessgeneral.h"
#include "lensgaconvergenceparameters.h"
#include "randomnumbergenerator.h"
#include <serut/memoryserializer.h>
#include <eatk/evolutionaryalgorithm.h>
#include <eatk/singlethreadedpopulationfitnesscalculation.h>
#include <eatk/multithreadedpopulationfitnesscalculation.h>
#include <eatk/multipopulationevolver.h>
#include <eatk/fasternondominatedsetcreator.h>
#include <eatk/fitnessbasedduplicateremoval.h>
#include <serut/vectorserializer.h>

#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <limits>
#include <chrono>

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
	// 	cout << "# Generation " << generation << ", calculated: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
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

class PreferredIndividualSelector
{
public:
	PreferredIndividualSelector() { }
	virtual ~PreferredIndividualSelector() { }
	virtual errut::bool_t select(const std::vector<std::shared_ptr<eatk::Individual>> &best,
								 std::shared_ptr<eatk::Individual> &selected)
	{
		return "Not implemented";
	}
};

class SubsequentBestIndividualSelector : public PreferredIndividualSelector
{
public:
	SubsequentBestIndividualSelector(size_t numObjectives,
									 const std::shared_ptr<eatk::FitnessComparison> &fitComp)
		: m_numObjectives(numObjectives), m_cmp(fitComp) { }
	~SubsequentBestIndividualSelector() { }
	
	errut::bool_t select(const std::vector<std::shared_ptr<eatk::Individual>> &bestIndividuals,
								 std::shared_ptr<eatk::Individual> &selected) override
	{
		std::vector<std::shared_ptr<eatk::Individual>> genomes = bestIndividuals;
		std::vector<std::shared_ptr<eatk::Individual>> genomes2;
		std::shared_ptr<eatk::Individual> best;

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
	std::shared_ptr<eatk::FitnessComparison> m_cmp;
};

class Stop : public grale::LensGAStopCriterion
{
public:
	Stop(const std::shared_ptr<grale::LensGAGenomeMutation> &mutation)
		: grale::LensGAStopCriterion(mutation) { }
protected:
	void onReport(const std::string &s)	const override
	{
		WriteLineStdout("GAMESSAGESTR:" + s);
	}
};

class MyExchange : public eatk::SequentialRandomIndividualExchange
{
public:
	MyExchange(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, size_t iterations) : eatk::SequentialRandomIndividualExchange(rng, iterations) { }
protected:
	void onExchange(size_t generation, size_t srcPop, size_t srcIndividualIdx, size_t dstPop, size_t dstIndividualIdx) override
	{
		std::cerr << "Generation " << generation << ": migrating " << srcIndividualIdx << " from pop " << srcPop << " to " << dstIndividualIdx << " in pop " << dstPop << std::endl;
	}
};

// TODO: rename this, is from copy-paste
class NewGACommunicatorBase : public InversionCommunicator
{
public:
	NewGACommunicatorBase() { }
	~NewGACommunicatorBase() { }

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
protected:	
	bool hasGeneticAlgorithm() const override { return true; }
	void getAllBestGenomes(std::vector<std::shared_ptr<eatk::Individual>> &bestGenomes) override
	{
		bestGenomes = m_best;
	}

	std::shared_ptr<eatk::Individual> getPreferredBestGenome() override
	{
		if (m_best.size() == 0)
			return nullptr;
		if (!m_selector.get())
			return nullptr;
		
		errut::bool_t r;
		std::shared_ptr<eatk::Individual> selected;
		if (!(r = m_selector->select(m_best, selected)))
		{
			std::cerr << "Error in preferred individual selection: " << r.getErrorString() << std::endl;
			return nullptr;
		}

		return selected;
	}

	errut::bool_t runGA(int popSize, const std::string &lensFitnessObjectType,
						 const std::string &calculatorType,
	                     grale::LensGACalculatorFactory &calcFactory, 
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
						 const std::vector<uint8_t> &factoryParamBytes,
						 const grale::GAParameters &params,
						 const grale::LensGAConvergenceParameters &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams)
	{
		errut::bool_t r;

		auto comparison = std::make_shared<grale::LensGAFitnessComparison>();

		m_selector = std::make_shared<SubsequentBestIndividualSelector>(
								genomeCalculator->getNumberOfObjectives(),
								comparison);

		std::shared_ptr<grale::RandomNumberGenerator> rng = std::make_shared<grale::RandomNumberGenerator>();
		MyGA ga;

		WriteLineStdout("GAMESSAGESTR:RNG SEED: " + std::to_string(rng->getSeed()));

		grale::LensGAIndividualCreation creation(rng, 
						  genomeCalculator->getNumberOfBasisFunctions(),
						  genomeCalculator->getNumberOfSheets(),
						  genomeCalculator->allowNegativeValues(),
						  genomeCalculator->getNumberOfObjectives());

		// TODO? For now, only the fitness calculations are parallel, using
		//       The same mutation instance for all populations is safe.
		//       In case a multi-threaded approach is used, this should be
		//       investigated again.
		//       The stop criterion is coupled to this (to switch from large
		//       mutations to small mutations for example), so we can't simply
		//       use multiple instances.
		auto mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 
					   1.0, // chance multiplier; has always been set to one
					   genomeCalculator->allowNegativeValues(),
					   0, // mutation amplitude, will be in the stop criterion
					   true); // absolute or small mutation, will be set in the stop criterion

		auto getEvolver = [&rng, &genomeCalculator, &params, &mutation]()
		{
			std::shared_ptr<grale::LensGACrossoverBase> cross;
			if (genomeCalculator->getNumberOfObjectives() == 1)
				cross = std::make_shared<grale::LensGASingleObjectiveCrossover>(params.getSelectionPressure(),
							  params.getUseElitism(),
							  params.getAlwaysIncludeBest(),
							  params.getCrossOverRate(),
							  rng,
							  genomeCalculator->allowNegativeValues(),
							  mutation);
			else
				cross = std::make_shared<grale::LensGAMultiObjectiveCrossover>(params.getSelectionPressure(),
							  params.getUseElitism(),
							  params.getAlwaysIncludeBest(),
							  params.getCrossOverRate(),
							  rng,
							  genomeCalculator->allowNegativeValues(),
							  mutation,
							  genomeCalculator->getNumberOfObjectives());

			return cross;
		};

		std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
		if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
								factoryParamBytes, creation, calc)))
			return "Can't get calculator: " + r.getErrorString();
		
		// Note: this approach causes the settings (large/small mutations) to switch
		// at the same time for all subpopulations. Is the alternative better?
		Stop stop(mutation);

		if (!(r = stop.initialize(genomeCalculator->getNumberOfObjectives(), convParams)))
		{
			calculatorCleanup();
			return "Error initializing convergence checker: " + r.getErrorString();
		}

		if (!multiPopParams.get()) // Single population only
		{
			auto cross = getEvolver();

			if (!(r = ga.run(creation, *cross, *calc, stop, popSize)))
			{
				calculatorCleanup();
				return "Error running GA: " + r.getErrorString();
			}

			m_best = cross->getBestIndividuals();
		}
		else // Use several populations, with migration
		{
			size_t numPop = multiPopParams->getNumberOfPopulations();
			if (numPop < 2)
			{
				calculatorCleanup();
				return "At least 2 populations are needed for a multi-population GA";
			}
			if (numPop > 64) // TODO: what's a reasonable upper limit?
			{
				calculatorCleanup();
				return "Currently there's a maximum of 64 populations";
			}

			std::vector<size_t> popSizes;
			for (size_t i = 0 ; i < numPop ; i++)
				popSizes.push_back(popSize);

			std::shared_ptr<eatk::BestIndividualMerger> merger;
			size_t numObjectives = genomeCalculator->getNumberOfObjectives();
			if (numObjectives == 1)
				merger = std::make_shared<eatk::SingleObjectiveBestIndividualMerger>(comparison);
			else
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				auto dupRemoval = std::make_shared<eatk::FitnessBasedDuplicateRemoval>(comparison, numObjectives);
				merger = std::make_shared<eatk::MultiObjectiveBestIndividualMerger>(ndCreator, dupRemoval);
			}

			std::vector<std::shared_ptr<eatk::PopulationEvolver>> evolvers;
			for (size_t i = 0 ; i < numPop ; i++)
				evolvers.push_back(getEvolver());

			eatk::MultiPopulationEvolver multiPopEvolver(evolvers, merger);

			auto migrationCheck = std::make_shared<eatk::UniformProbabilityMigrationCheck>(rng,
					                                                                       (float)multiPopParams->getMigrationGenerationFraction(),
																						   multiPopParams->getNumberOfInitialPopulationsToSkip());
			auto migrationExchange = std::make_shared<MyExchange>(rng, multiPopParams->getNumberOfIndividualsToLeavePopulation());

			eatk::BasicMigrationStrategy migration(migrationCheck, migrationExchange);

			if (!(r = ga.run(creation, multiPopEvolver, *calc, stop, migration, popSizes)))
			{
				calculatorCleanup();
				return "Error running GA: " + r.getErrorString();
			}

			m_best = multiPopEvolver.getBestIndividuals();
		}

		calculatorCleanup();

		// std::cout << "Best: " << std::endl;
		// for (auto &b: m_best)
		// 	std::cout << b->fitness()->toString() << std::endl;

		return true;
	}

	virtual errut::bool_t getCalculator(const std::string &lensFitnessObjectType,
									const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									grale::LensGAIndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) = 0;

	virtual void calculatorCleanup() { }
	
	std::vector<std::shared_ptr<eatk::Individual>> m_best;
	std::shared_ptr<SubsequentBestIndividualSelector> m_selector;
};
