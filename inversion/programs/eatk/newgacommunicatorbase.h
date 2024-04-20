#pragma once

#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "deparameters.h"
#include "lensinversiongafactorycommon.h"
#include "inversioncommunicatornewga.h"
#include "constants.h"
#include "lensgaindividual.h"
#include "lensgagenomecrossover.h"
#include "lensgafitnesscomparison.h"
#include "lensgasingleobjectivecrossover.h"
#include "lensgamultiobjectivecrossover.h"
#include "lensdemutationcrossover.h"
#include "lensdeevolver.h"
#include "lensfitnessgeneral.h"
#include "lensgaconvergenceparameters.h"
#include "randomnumbergenerator.h"
#include "ea.h"
#include "preferredindividualselector.h"
#include "stop.h"
#include "exchange.h"
#include "reusecreation.h"
#include <serut/memoryserializer.h>
#include <eatk/singlethreadedpopulationfitnesscalculation.h>
#include <eatk/multithreadedpopulationfitnesscalculation.h>
#include <eatk/fasternondominatedsetcreator.h>
#include <eatk/fitnessbasedduplicateremoval.h>
#include <serut/vectorserializer.h>

#include <sstream>
#include <string>
#include <limits>

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

		if (eaType == "GA" || eaType == "GA+JADE")
			r = runGA_GA(rng, creation, comparison, *calc, popSize, genomeCalculator->allowNegativeValues(), genomeCalculator->getNumberOfObjectives(),
					        params, convParams, multiPopParams, eaType);
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
						 const std::string &eaType, size_t generationOffsetForReporting = 0)
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

		Stop stop(mutation, -1, generationOffsetForReporting);
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
			bool needStrictlyBetter = params.getNeedStrictlyBetter();

			WriteLineStdout("GAMESSAGESTR:Running JADE algorithm, p = " + std::to_string(p) + 
					        ", c = " + std::to_string(c) + ", useArchive = " + std::to_string((int)useArch) +
							", initMuF = " + std::to_string(initMuF) + ", initMuCR = " + std::to_string(initMuCR) + 
							", needStrictlyBetter = " + std::to_string(needStrictlyBetter));

			if (numObjectives == 1)
			{
				evolver = std::make_unique<grale::LensJADEEvolver>(rng, mut, cross, comparison, 0,
						                                           p, c, useArch, initMuF, initMuCR,
																   1, nullptr, needStrictlyBetter);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<grale::LensJADEEvolver>(rng, mut, cross, comparison,
						  -1, // signals multi-objective
						  p, c, useArch, initMuF, initMuCR,
						  numObjectives, ndCreator,
						  needStrictlyBetter);
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
			bool needStrictlyBetter = params.getNeedStrictlyBetter();

			WriteLineStdout("GAMESSAGESTR:Running DE algorithm, F = " + std::to_string(F) + ", CR = " + std::to_string(CR) +
					        ", needStrictlyBetter = " + std::to_string(needStrictlyBetter));

			if (numObjectives == 1) // Single objective
			{
				evolver = std::make_unique<grale::LensDEEvolver>(rng, mut, F, cross, CR, comparison,
						                                         0, 1, nullptr, needStrictlyBetter);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<grale::LensDEEvolver>(rng, mut, params.getF(), cross, params.getCR(), comparison,
						                                         -1, numObjectives, ndCreator,
																 needStrictlyBetter); // -1 signals multi-objective
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
						 const grale::LensGAConvergenceParameters &convParams0,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::string &eaType)
	{
		const grale::GAParameters *pParams = dynamic_cast<const grale::GAParameters*>(&eaParams);
		if (!pParams)
			return "Invalid EA parameters for GA";
		const grale::GAParameters &params = *pParams;

		WriteLineStdout("GAMESSAGESTR:Running GA algorithm, selection pressure = " + std::to_string(params.getSelectionPressure()) +
				        ", elitism = " + std::to_string((int)params.getUseElitism()) +
						", always include best = " + std::to_string((int)params.getAlwaysIncludeBest()) + 
						", crossover rate = " + std::to_string(params.getCrossOverRate()));


		grale::LensGAConvergenceParameters convParamsGA;
		grale::LensGAConvergenceParameters convParamsJADE;

		if (eaType == "GA")
		{
			// Nothing to do
			convParamsGA = convParams0;
		}
		else if (eaType == "GA+JADE")
		{
			WriteLineStdout("GAMESSAGESTR:Will add JADE after standard GA run");
			if (multiPopParams.get())
				return "GA+JADE doesn't work with multiple populations";

			// To make sure that history size and max generations are copies
			convParamsGA = convParams0;
			convParamsJADE = convParams0;

			// Split conversion parameters: last is for JADE
			std::vector<double> factors = convParams0.getConvergenceFactors();
			std::vector<double> mutSizes = convParams0.getConvergenceSmallMutationSizes();

			if (factors.size() > 0 && mutSizes.size() > 0)
			{
				double deConvFactor = factors.back();
				double deMutSize = mutSizes.back(); // This isn't actually used though
				factors.resize(factors.size()-1); // Remove last
				mutSizes.resize(mutSizes.size()-1);

				if (!convParamsGA.setConvergenceFactorsAndMutationSizes(factors, mutSizes))
					return "Can't set GA convergence factors: " + convParamsGA.getErrorString();
				if (!convParamsJADE.setConvergenceFactorsAndMutationSizes({deConvFactor}, {deMutSize}))
					return "Can't set JADE convergence factors: " + convParamsJADE.getErrorString();
			}
		}
		else
			return "Unexpected eaType '" + eaType + "'";

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

			if (!(r = stop.initialize(numObjectives, convParamsGA)))
				return "Error initializing convergence checker: " + r.getErrorString();

			if (!(r = ga.run(creation, *cross, calc, stop, popSize)))
				return "Error running GA: " + r.getErrorString();

			if (eaType == "GA+JADE")
			{
				WriteLineStdout("GAMESSAGESTR:EXPERIMENTAL JADE FINISH");

				ReuseCreation reuseCreation(*(ga.getPopulation()));
				size_t generationOffset = ga.getNumberOfGenerations();

				grale::JADEParameters defaultJADEParams; // TODO: make this configurable somehow

				if (!(r = runGA_DE(rng, reuseCreation, comparison, calc, popSize, allowNegative, numObjectives, defaultJADEParams,
								   convParamsJADE, nullptr, "JADE", generationOffset)))
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
			if (!(r = stop.initialize(numObjectives, convParamsGA)))
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
