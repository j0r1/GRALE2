#pragma once

#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "deparameters.h"
#include "rndparameters.h"
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
#include "rngwrapper.h"
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
						 const std::vector<std::unique_ptr<grale::EAParameters>> &allEAParams,
						 const std::vector<grale::LensGAConvergenceParameters> &allConvParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::vector<std::string> &allEATypes)
	{
		if (allEATypes.size() != allEAParams.size() || allEATypes.size() != allConvParams.size())
			return "Unexpected mismatch between number of EA types (" + std::to_string(allEATypes.size()) +
				   "), EA parameters (" + std::to_string(allEAParams.size()) + ") and convergence parameters (" + 
				   std::to_string(allConvParams.size()) + ")";

		// TODO: Check type name and parameters
		for (size_t i = 0 ; i < allEATypes.size() ; i++)
		{
			// TODO
		}


		std::shared_ptr<grale::RandomNumberGenerator> rng0 = std::make_shared<grale::RandomNumberGenerator>();
		WriteLineStdout("GAMESSAGESTR:RNG SEED: " + std::to_string(rng0->getSeed()));

#if 0
		std::shared_ptr<eatk::RandomNumberGenerator> rng = std::make_shared<RngWrapper>(rng0);
#else
		std::shared_ptr<eatk::RandomNumberGenerator> rng = rng0;
#endif


		auto comparison = std::make_shared<grale::LensGAFitnessComparison>();
		m_selector = std::make_shared<SubsequentBestIndividualSelector>(
								genomeCalculator->getNumberOfObjectives(),
								comparison);

		// After this is created, calculatorCleanup() should be called as well (for MPI at the moment)
		std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
		errut::bool_t r;

		std::unique_ptr<eatk::IndividualCreation> creation;
		{
			std::unique_ptr<grale::LensGAIndividualCreation> lensGACreation = std::make_unique<grale::LensGAIndividualCreation>(rng, 
							  genomeCalculator->getNumberOfBasisFunctions(),
							  genomeCalculator->getNumberOfSheets(),
							  genomeCalculator->allowNegativeValues(),
							  genomeCalculator->getNumberOfObjectives());

			if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
									factoryParamBytes, *lensGACreation, calc)))
				return "Can't get calculator: " + r.getErrorString();

			creation = std::move(lensGACreation);
		}

		// For compatibility with previous approach, in multi-objective GA we need to copy this as
		// well (for elitism); it's a vector of vectors for the multi-pop case
		std::vector<std::vector<std::shared_ptr<eatk::Individual>>> previousBestIndividuals;
		size_t generationCount = 0;

		for (size_t i = 0 ; i < allEATypes.size() ; i++)
		{
			// Don't process this anymore
			if (generationCount > allConvParams[i].getMaximumNumberOfGenerations())
			{
				WriteLineStdout("GAMESSAGESTR:DEBUG: generationCount (" + std::to_string(generationCount) + ") > allConvParams[" + std::to_string(i) + 
						        "].getMaximumNumberOfGenerations() (" + std::to_string(allConvParams[i].getMaximumNumberOfGenerations()) + "), skipping next EA in line (" + allEATypes[i] + ")");
				continue;
			}

			std::string eaType = allEATypes[i];
			size_t numGen = 0;
			std::unique_ptr<eatk::IndividualCreation> reuseCreation;

			if (eaType == "GA")
				std::tie(r, previousBestIndividuals, reuseCreation, numGen) = runGA_GA(rng, *creation, comparison, *calc, popSize, genomeCalculator->allowNegativeValues(), genomeCalculator->getNumberOfObjectives(),
								*(allEAParams[i]), allConvParams[i], multiPopParams, eaType, generationCount, previousBestIndividuals);
			else if (eaType == "DE" || eaType == "JADE")
				std::tie(r, previousBestIndividuals, reuseCreation, numGen) = runGA_DE(rng, *creation, comparison, *calc, popSize, genomeCalculator->allowNegativeValues(), genomeCalculator->getNumberOfObjectives(),
								*(allEAParams[i]), allConvParams[i], multiPopParams, eaType, generationCount, previousBestIndividuals);
			else if (eaType == "RND")
				std::tie(r, previousBestIndividuals, reuseCreation, numGen) = runGA_RND(rng, *creation, popSize, 
								*(allEAParams[i]), multiPopParams);
			else
				r = "Unknown EA type '" + eaType + "', should be either GA, DE or JADE";

			if (!r)
				break;

			std::swap(creation, reuseCreation); // Make sure the next algorithm starts where this one left off

			generationCount += numGen;
		}

		// Note: m_best must be set inside the subroutines; previousBest can be the one from multiple
		//       populations, don't want to recalculate the non-dominated set here

		calculatorCleanup();
		return r;
	}

	static inline std::tuple<errut::bool_t,
		                                std::vector<std::vector<std::shared_ptr<eatk::Individual>>>,
						                std::unique_ptr<eatk::IndividualCreation>, size_t> E(const std::string &msg)

	{
		std::vector<std::vector<std::shared_ptr<eatk::Individual>>> dummy;
		errut::bool_t r = msg;
		return { r, dummy, nullptr, 0 };
	};

	std::tuple<errut::bool_t,
		       std::vector<std::vector<std::shared_ptr<eatk::Individual>>>,
			   std::unique_ptr<eatk::IndividualCreation>,
			   size_t> runGA_RND(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 eatk::IndividualCreation &creation,
			             int popSize,
						 const grale::EAParameters &eaParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams)
	{
		if (dynamic_cast<ReuseCreation*>(&creation) == nullptr)
			return E("Seems that RND step is first in line, but it should use the result from another algorithm");

		const grale::RNDParameters *pParams = dynamic_cast<const grale::RNDParameters*>(&eaParams);
		if (pParams == nullptr)
			return E("Expecting parameters of type RNDParameters for RND step");
		const grale::RNDParameters &params = *pParams;
		double scale = params.getScale();
		double minFactor = 1.0-scale;
		double maxFactor = 1.0+scale;

		if (minFactor < 0)
			return E("Got negative factor from scale factor " + std::to_string(minFactor));

		WriteLineStdout("GAMESSAGESTR:Running RND, applying some random factor between 1-scale and 1+scale, scale = " + std::to_string(scale));
		WriteLineStdout("GAMESSAGESTR:Note that this RND step will not update the best solutions so far");

		// We'll use a temporary population to store the new genomes in
		std::shared_ptr<eatk::Population> tmpPop = std::make_shared<eatk::Population>();
		std::shared_ptr<eatk::Genome> genome;
		while (true)
		{
			genome = creation.createInitializedGenome();
			if (genome.get() == nullptr) // We've processed everything
				break;

			// Create copy and modify it slightly
			genome = genome->createCopy();
			grale::LensGAGenome &g = static_cast<grale::LensGAGenome&>(*genome);
			g.m_parent1 = -1;
			g.m_parent2 = -1;

			for (auto &x : g.m_weights)
				x *= rng->getRandomFloat((float)minFactor, (float)maxFactor);
			for (auto &x : g.m_sheets)
				x *= rng->getRandomFloat((float)minFactor, (float)maxFactor);

			// Add this to our new population
			auto ind = std::make_shared<grale::LensGAIndividual>(genome, creation.createEmptyFitness());
			tmpPop->append(ind);
		}

		// Note that we're not changing m_best in this step!

		std::vector<std::shared_ptr<eatk::Population>> populations { tmpPop };
		std::vector<std::vector<std::shared_ptr<eatk::Individual>>> emptyBest;

		return { true, emptyBest, std::make_unique<ReuseCreation>(populations), 0 };

	}

	std::tuple<errut::bool_t,
		       std::vector<std::vector<std::shared_ptr<eatk::Individual>>>,
			   std::unique_ptr<eatk::IndividualCreation>,
			   size_t> runGA_DE(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 eatk::IndividualCreation &creation,
						 const std::shared_ptr<grale::LensGAFitnessComparison> &comparison,
						 eatk::PopulationFitnessCalculation &calc,
			             int popSize,
						 bool allowNegative, size_t numObjectives,
						 const grale::EAParameters &eaParams,
						 const grale::LensGAConvergenceParameters &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::string &eaType, size_t generationOffsetForReporting,
						 const std::vector<std::vector<std::shared_ptr<eatk::Individual>>> &previousBest)
	{
		if (previousBest.size() > 0)
			WriteLineStdout("GAMESSAGESTR:WARNING: nog using previous best for initialization of DE/JADE");

		errut::bool_t r;

		if (multiPopParams.get())
			return E("DE/JADE only works with a single population");

		Stop stop(-1, generationOffsetForReporting);
		if (!(r = stop.initialize(numObjectives, convParams)))
			return E("Can't initialize stop criterion: " + r.getErrorString());

		auto mut = std::make_shared<grale::LensDEMutation>();
		auto cross = std::make_shared<grale::LensDECrossover>(rng, allowNegative);

		std::unique_ptr<eatk::PopulationEvolver> evolver;

		if (eaType == "JADE")
		{
			const grale::JADEParameters *pParams = dynamic_cast<const grale::JADEParameters*>(&eaParams);
			if (!pParams)
				return E("Invalid EA parameters for JADE");

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
				return E("Invalid EA parameters for DE");

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
		else
			return E("Unexpected eaType '" + eaType + "'");

		MyGA ga;
		if (!(r = ga.run(creation, *evolver, calc, stop, popSize, popSize, popSize*2)))
			return E("Error running GA: " + r.getErrorString());

		std::vector<std::vector<std::shared_ptr<eatk::Individual>>> allBest = { { evolver->getBestIndividuals() } };
		m_best = evolver->getBestIndividuals();

		return { true, allBest,
				 std::make_unique<ReuseCreation>(ga.getPopulations()),
				 ga.getNumberOfGenerations() };

	}

	std::tuple<errut::bool_t,
		       std::vector<std::vector<std::shared_ptr<eatk::Individual>>>,
			   std::unique_ptr<eatk::IndividualCreation>,
			   size_t> runGA_GA(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 eatk::IndividualCreation &creation,
						 const std::shared_ptr<grale::LensGAFitnessComparison> &comparison,
						 eatk::PopulationFitnessCalculation &calc,
			             int popSize,
						 bool allowNegative, size_t numObjectives,
						 const grale::EAParameters &eaParams,
						 const grale::LensGAConvergenceParameters &convParamsGA,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::string &eaType, size_t generationOffsetForReporting,
						 const std::vector<std::vector<std::shared_ptr<eatk::Individual>>> &previousBestIndividuals)
	{
		const grale::GAParameters *pParams = dynamic_cast<const grale::GAParameters*>(&eaParams);
		if (!pParams)
			return E("Invalid EA parameters for GA");
		const grale::GAParameters &params = *pParams;

		WriteLineStdout("GAMESSAGESTR:Running GA algorithm, selection pressure = " + std::to_string(params.getSelectionPressure()) +
				        ", elitism = " + std::to_string((int)params.getUseElitism()) +
						", always include best = " + std::to_string((int)params.getAlwaysIncludeBest()) + 
						", crossover rate = " + std::to_string(params.getCrossOverRate()) +
						", small mutation size = " + std::to_string(params.getSmallMutationSize()));


		errut::bool_t r;
		MyGA ga;

		double smallMutSize = params.getSmallMutationSize();
		bool absoluteMutation = (smallMutSize <= 0)?true:false;

		auto mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 
					   1.0, // chance multiplier; has always been set to one
					   allowNegative,
					   smallMutSize,
					   absoluteMutation);

		auto getEvolver = [&rng, allowNegative, numObjectives, &params, &mutation]()
		{
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

		auto restorePreviousBest = [&previousBestIndividuals](std::shared_ptr<grale::LensGACrossoverBase> &evolver, size_t i)
		{
			assert(i < previousBestIndividuals.size());

			if (auto moCross = dynamic_cast<grale::LensGAMultiObjectiveCrossover*>(evolver.get()) ; moCross != nullptr)
			{
				WriteLineStdout("GAMESSAGESTR:DEBUG: multi-objective GA, restoring best " + std::to_string(previousBestIndividuals[i].size()) + " individuals");
				moCross->restoreBestIndividuals(previousBestIndividuals[i]);
			}
			else
				WriteLineStdout("GAMESSAGESTR:DEBUG: not a multi-objective GA, ignoring restoring best " + std::to_string(previousBestIndividuals[i].size()) + " individuals");
		};

		if (!multiPopParams.get()) // Single population only
		{
			auto cross = getEvolver();

			if (previousBestIndividuals.size() > 1) // expecting one population
				return E("Only expecting previous best individuals from previous EA step for one population, but got " + std::to_string(previousBestIndividuals.size()));

			if (previousBestIndividuals.size() == 1)
				restorePreviousBest(cross, 0);
			else
				WriteLineStdout("GAMESSAGESTR:DEBUG: no previous best to restore in GA");

			Stop stop(-1, generationOffsetForReporting);

			if (!(r = stop.initialize(numObjectives, convParamsGA)))
				return E("Error initializing convergence checker: " + r.getErrorString());

			if (!(r = ga.run(creation, *cross, calc, stop, popSize)))
				return E("Error running GA: " + r.getErrorString());

			m_best = cross->getBestIndividuals();

			std::vector<std::vector<std::shared_ptr<eatk::Individual>>> allBest = { { cross->getBestIndividuals() } };

			return { true, allBest,
				 std::make_unique<ReuseCreation>(ga.getPopulations()),
				 ga.getNumberOfGenerations() };

		}
		else // Use several populations, with migration
		{
			size_t numPop = multiPopParams->getNumberOfPopulations();
			if (numPop < 2)
				return E("At least 2 populations are needed for a multi-population GA");

			if (numPop > 64) // TODO: what's a reasonable upper limit?
				return E("Currently there's a maximum of 64 populations");

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

			if (previousBestIndividuals.size() == 0)
				WriteLineStdout("GAMESSAGESTR:DEBUG: no previous best to restore in multi-population GA");
			else if (previousBestIndividuals.size() != numPop)
				return E("Expecting to restore previous best individuals from " + std::to_string(numPop) + " subpopulations, but got " + std::to_string(previousBestIndividuals.size()));

			std::vector<std::shared_ptr<eatk::PopulationEvolver>> evolvers;
			for (size_t i = 0 ; i < numPop ; i++)
			{
				auto evolver = getEvolver();
				evolvers.push_back(evolver);

				// Restore previous best for each evolver
				if (previousBestIndividuals.size() > 0) // We've already checked the size matches numPop
					restorePreviousBest(evolver, i);
			}

			eatk::MultiPopulationEvolver multiPopEvolver(evolvers, merger);

			// Adjust numberOfInitialGenerationsToSkip based on generationOffsetForReporting
			size_t numGenerationsToSkip = multiPopParams->getNumberOfInitialGenerationsToSkip();
			if (generationOffsetForReporting > numGenerationsToSkip)
				numGenerationsToSkip = 0;
			else
				numGenerationsToSkip -= generationOffsetForReporting;

			if (numGenerationsToSkip != multiPopParams->getNumberOfInitialGenerationsToSkip())
				WriteLineStdout("GAMESSAGESTR:DEBUG: adjusted numGenerationsToSkip from " + 
						        std::to_string(multiPopParams->getNumberOfInitialGenerationsToSkip()) + " to " + 
								std::to_string(numGenerationsToSkip) + " based on generationOffsetForReporting (" + 
								std::to_string(generationOffsetForReporting) + ")");

			auto migrationCheck = std::make_shared<eatk::UniformProbabilityMigrationCheck>(rng,
					                                                                       (float)multiPopParams->getMigrationGenerationFraction(),
																						   numGenerationsToSkip);
			auto migrationExchange = std::make_shared<MyExchange>(rng, multiPopParams->getNumberOfIndividualsToLeavePopulation());

			eatk::BasicMigrationStrategy migration(migrationCheck, migrationExchange);

			MultiStop stop(evolvers.size(), generationOffsetForReporting);
			if (!(r = stop.initialize(numObjectives, convParamsGA)))
				return E("Error initializing multi-population convergence checker: " + r.getErrorString());

			if (!(r = ga.run(creation, multiPopEvolver, calc, stop, migration, popSizes)))
				return E("Error running GA: " + r.getErrorString());

			m_best = multiPopEvolver.getBestIndividuals();

			std::vector<std::vector<std::shared_ptr<eatk::Individual>>> allBest;
			for (auto &evolver : evolvers)
				allBest.push_back(evolver->getBestIndividuals());

			return { true, allBest,
				 std::make_unique<ReuseCreation>(ga.getPopulations()),
				 ga.getNumberOfGenerations() };
		}
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
