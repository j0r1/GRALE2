#pragma once

#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "deparameters.h"
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
#include "lensdemutationcrossover.h"
#include "lensdeevolver.h"
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

class MultiStop : public eatk::StopCriterion
{
public:
	MultiStop(const std::vector<std::shared_ptr<grale::LensGAGenomeMutation>> &mutations)
	{
		for (int i = 0 ; i < (int)mutations.size() ; i++)
			m_stops.push_back(std::make_shared<Stop>(mutations[i], i));
	}

	errut::bool_t initialize(size_t numObjectives, const grale::LensGAConvergenceParameters &convParams)
	{
		errut::bool_t r;
		for (auto &stop : m_stops)
		{
			if (!(r = stop->initialize(numObjectives, convParams)))
				return "Unable to initialize stop criterion for subpopulation: " + r.getErrorString();
		}
		return true;
	}

	errut::bool_t analyze(const eatk::PopulationEvolver &ev, size_t generation, bool &shouldStop)
	{
		if (generation <= 1)
		{
			if (dynamic_cast<const eatk::MultiPopulationEvolver *>(&ev) == nullptr)
				return "Evolver doesn't appear to be a 'MultiPopulationEvolver'";
		}

		const eatk::MultiPopulationEvolver &multiEvolver = static_cast<const eatk::MultiPopulationEvolver &>(ev);
		auto &singleEvolvers = multiEvolver.getSinglePopulationEvolvers();

		if (singleEvolvers.size() != m_stops.size())
			return "Number of single population evolvers (" + std::to_string(singleEvolvers.size()) + ") doesn't match number individual stop criteria (" + std::to_string(m_stops.size()) + ")";

		bool stop = true;
		errut::bool_t r;

		for (size_t i = 0 ; i < m_stops.size() ; i++)
		{
			bool shouldStopSingle = false;

			if (!(r = m_stops[i]->analyze(*singleEvolvers[i], generation, shouldStopSingle)))
				return "Error running stop criterion for population " + std::to_string(i) + ": " + r.getErrorString();
			if (!shouldStopSingle)
				stop = false;
		}

		// Stop if all populations indicate stop
		shouldStop = stop;

		return true;
	}
private:
	std::vector<std::shared_ptr<Stop>> m_stops;
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

class ReuseCreation : public eatk::IndividualCreation
{
public:
	ReuseCreation(const eatk::Population &pop)
	{
		if (pop.size() == 0)
			return;

		m_referenceIndividual = pop.individual(0)->createCopy();
		m_referenceGenome = m_referenceIndividual->genomePtr()->createCopy();
		m_referenceFitness = m_referenceIndividual->fitnessPtr()->createCopy();

		for (auto &i : pop.individuals())
			m_genomePool.push_back(i->genomePtr()->createCopy());
	}

	std::shared_ptr<eatk::Genome> createInitializedGenome() override
	{
		if (m_genomePool.size() == 0)
			return nullptr;
		auto genome = m_genomePool.front();
		m_genomePool.pop_front();
		return genome;
	}

	std::shared_ptr<eatk::Genome> createUnInitializedGenome() override
	{ 
		if (!m_referenceGenome.get())
			return nullptr;
		return m_referenceGenome->createCopy(false);
	}

	std::shared_ptr<eatk::Fitness> createEmptyFitness() override
	{
		if (!m_referenceFitness.get())
			return nullptr;
		return m_referenceFitness->createCopy(false);
	}

	std::shared_ptr<eatk::Individual> createReferenceIndividual() override { return m_referenceIndividual; }
private:
	std::list<std::shared_ptr<eatk::Genome>> m_genomePool;
	std::shared_ptr<eatk::Individual> m_referenceIndividual;
	std::shared_ptr<eatk::Genome> m_referenceGenome;
	std::shared_ptr<eatk::Fitness> m_referenceFitness;
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
						 const grale::EAParameters &params,
						 const grale::LensGAConvergenceParameters &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::string &eaType)
	{
		std::shared_ptr<grale::RandomNumberGenerator> rng = std::make_shared<grale::RandomNumberGenerator>();
		WriteLineStdout("GAMESSAGESTR:RNG SEED: " + std::to_string(rng->getSeed()));

		grale::LensGAIndividualCreation creation(rng, 
						  genomeCalculator->getNumberOfBasisFunctions(),
						  genomeCalculator->getNumberOfSheets(),
						  genomeCalculator->allowNegativeValues(),
						  genomeCalculator->getNumberOfObjectives());

		auto comparison = std::make_shared<grale::LensGAFitnessComparison>();
		m_selector = std::make_shared<SubsequentBestIndividualSelector>(
								genomeCalculator->getNumberOfObjectives(),
								comparison);

		// After this is created, calculatorCleanup() should be called as well (for MPI at the moment)
		std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
		errut::bool_t r;

		if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
								factoryParamBytes, creation, calc)))
			return "Can't get calculator: " + r.getErrorString();

		if (eaType == "GA")
			r = runGA_GA(rng, creation, comparison, *calc, popSize, genomeCalculator->allowNegativeValues(), genomeCalculator->getNumberOfObjectives(),
					        params, convParams, multiPopParams);
		else if (eaType == "DE" || eaType == "JADE")
			r = runGA_DE(rng, creation, comparison, *calc, popSize, genomeCalculator->allowNegativeValues(), genomeCalculator->getNumberOfObjectives(),
					        params, convParams, multiPopParams, eaType);
		else
			r = "Unknown EA type '" + eaType + "', should be either GA, DE or JADE";

		calculatorCleanup();
		return r;
	}

	// TODO: merge more common code from the GA and DE versions

	errut::bool_t runGA_DE(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 eatk::IndividualCreation &creation,
						 const std::shared_ptr<grale::LensGAFitnessComparison> &comparison,
						 eatk::PopulationFitnessCalculation &calc,
			             int popSize,
						 bool allowNegative, size_t numObjectives,
						 const grale::EAParameters &eaParams,
						 const grale::LensGAConvergenceParameters &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::string &eaType)
	{
		errut::bool_t r;

		if (multiPopParams.get())
			return "DE/JADE only works with a single population";

		// TODO: At the moment we need this for the stop criterion
		std::shared_ptr<grale::LensGAGenomeMutation> mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 
						   1.0, // chance multiplier; has always been set to one
						   allowNegative,
						   0, // mutation amplitude, will be in the stop criterion
						   true); // absolute or small mutation, will be set in the stop criterion

		// TODO: this stop criterion where the mutation amplitude is changed doesn't really
		//       make sense (not the kind of mutation in DE)

		Stop stop(mutation);
		if (!(r = stop.initialize(numObjectives, convParams)))
			return "Can't initialize stop criterion: " + r.getErrorString();

		auto mut = std::make_shared<grale::LensDEMutation>();
		auto cross = std::make_shared<grale::LensDECrossover>(rng, allowNegative);

		std::unique_ptr<eatk::PopulationEvolver> evolver;

		if (eaType == "JADE")
		{
			const grale::JADEParameters *pParams = dynamic_cast<const grale::JADEParameters*>(&eaParams);
			if (!pParams)
				return "Invalid EA parameters for JADE";
			const grale::JADEParameters &params = *pParams;
			double p = params.getBestFraction_p();
			double c = params.getParameterUpdateFraction_c();
			bool useArch = params.useExternalArchive();
			double initMuF = params.getInitialMeanF();
			double initMuCR = params.getInitialMeanCR();

			WriteLineStdout("GAMESSAGESTR:Running JADE algorithm, p = " + std::to_string(p) + 
					        ", c = " + std::to_string(c) + ", useArchive = " + std::to_string((int)useArch) +
							", initMuF = " + std::to_string(initMuF) + ", initMuCR = " + std::to_string(initMuCR));

			if (numObjectives == 1)
			{
				evolver = std::make_unique<grale::LensJADEEvolver>(rng, mut, cross, comparison, 0,
						                                           p, c, useArch, initMuF, initMuCR);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<grale::LensJADEEvolver>(rng, mut, cross, comparison,
						  -1, // signals multi-objective
						  p, c, useArch, initMuF, initMuCR,
						  numObjectives, ndCreator);
			}
		}
		else if (eaType == "DE")
		{
			const grale::DEParameters *pParams = dynamic_cast<const grale::DEParameters*>(&eaParams);
			if (!pParams)
				return "Invalid EA parameters for DE";
			const grale::DEParameters &params = *pParams;
			double F = params.getF();
			double CR = params.getCR();

			WriteLineStdout("GAMESSAGESTR:Running DE algorithm, F = " + std::to_string(F) + ", CR = " + std::to_string(CR));

			if (numObjectives == 1) // Single objective
			{
				evolver = std::make_unique<grale::LensDEEvolver>(rng, mut, F, cross, CR, comparison);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<grale::LensDEEvolver>(rng, mut, params.getF(), cross, params.getCR(), comparison,
						                                         -1, numObjectives, ndCreator); // -1 signals multi-objective
			}
		}

		MyGA ga;
		if (!(r = ga.run(creation, *evolver, calc, stop, popSize, popSize, popSize*2)))
			return "Error running GA: " + r.getErrorString();

		m_best = evolver->getBestIndividuals();

		return true;
	}

	errut::bool_t runGA_GA(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 eatk::IndividualCreation &creation,
						 const std::shared_ptr<grale::LensGAFitnessComparison> &comparison,
						 eatk::PopulationFitnessCalculation &calc,
			             int popSize,
						 bool allowNegative, size_t numObjectives,
						 const grale::EAParameters &eaParams,
						 const grale::LensGAConvergenceParameters &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams)
	{
		const grale::GAParameters *pParams = dynamic_cast<const grale::GAParameters*>(&eaParams);
		if (!pParams)
			return "Invalid EA parameters for GA";
		const grale::GAParameters &params = *pParams;

		WriteLineStdout("GAMESSAGESTR:Running GA algorithm, selection pressure = " + std::to_string(params.getSelectionPressure()) +
				        ", elitism = " + std::to_string((int)params.getUseElitism()) +
						", always include best = " + std::to_string((int)params.getAlwaysIncludeBest()) + 
						", crossover rate = " + std::to_string(params.getCrossOverRate()));

		errut::bool_t r;
		MyGA ga;

		std::vector<std::shared_ptr<grale::LensGAGenomeMutation>> mutations;

		auto getEvolver = [&rng, allowNegative, numObjectives, &params, &mutations]()
		{
			auto mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 
						   1.0, // chance multiplier; has always been set to one
						   allowNegative,
						   0, // mutation amplitude, will be in the stop criterion
						   true); // absolute or small mutation, will be set in the stop criterion

			mutations.push_back(mutation);

			std::shared_ptr<grale::LensGACrossoverBase> cross;
			if (numObjectives == 1)
				cross = std::make_shared<grale::LensGASingleObjectiveCrossover>(params.getSelectionPressure(),
							  params.getUseElitism(),
							  params.getAlwaysIncludeBest(),
							  params.getCrossOverRate(),
							  rng,
							  allowNegative,
							  mutation);
			else
				cross = std::make_shared<grale::LensGAMultiObjectiveCrossover>(params.getSelectionPressure(),
							  params.getUseElitism(),
							  params.getAlwaysIncludeBest(),
							  params.getCrossOverRate(),
							  rng,
							  allowNegative,
							  mutation,
							  numObjectives);

			return cross;
		};

		if (!multiPopParams.get()) // Single population only
		{
			auto cross = getEvolver();

			assert(mutations.size() == 1);
			Stop stop(mutations[0]);

			if (!(r = stop.initialize(numObjectives, convParams)))
				return "Error initializing convergence checker: " + r.getErrorString();

			if (!(r = ga.run(creation, *cross, calc, stop, popSize)))
				return "Error running GA: " + r.getErrorString();

			if (std::getenv("GRALE_JADEFINISH"))
			{
				WriteLineStdout("GAMESSAGESTR:EXPERIMENTAL JADE FINISH");

				ReuseCreation reuseCreation(*(ga.getPopulation()));
				size_t generationOffset = ga.getNumberOfGenerations(); // TODO: use this in reporting

				grale::JADEParameters defaultJADEParams;
				grale::LensGAConvergenceParameters finalConvParams;

				finalConvParams.setConvergenceFactorsAndMutationSizes({0.05}, {-1});

				if (!(r = runGA_DE(rng, reuseCreation, comparison, calc, popSize, allowNegative, numObjectives, defaultJADEParams,
								   finalConvParams, nullptr, "JADE")))
					return "Can't run JADE finish: " + r.getErrorString();
			}
			else
				m_best = cross->getBestIndividuals();
		}
		else // Use several populations, with migration
		{
			size_t numPop = multiPopParams->getNumberOfPopulations();
			if (numPop < 2)
				return "At least 2 populations are needed for a multi-population GA";

			if (numPop > 64) // TODO: what's a reasonable upper limit?
				return "Currently there's a maximum of 64 populations";

			std::vector<size_t> popSizes;
			for (size_t i = 0 ; i < numPop ; i++)
				popSizes.push_back(popSize);

			std::shared_ptr<eatk::BestIndividualMerger> merger;
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
																						   multiPopParams->getNumberOfInitialGenerationsToSkip());
			auto migrationExchange = std::make_shared<MyExchange>(rng, multiPopParams->getNumberOfIndividualsToLeavePopulation());

			eatk::BasicMigrationStrategy migration(migrationCheck, migrationExchange);

			assert(mutations.size() > 1);
			MultiStop stop(mutations);
			if (!(r = stop.initialize(numObjectives, convParams)))
				return "Error initializing multi-population convergence checker: " + r.getErrorString();

			if (!(r = ga.run(creation, multiPopEvolver, calc, stop, migration, popSizes)))
				return "Error running GA: " + r.getErrorString();

			m_best = multiPopEvolver.getBestIndividuals();
		}

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
