#pragma once

#include "inversioncommunicatorbase.h"

class FreeFormInversionCommunicator : public InversionCommunicatorBase
{
public:
	FreeFormInversionCommunicator() { }
	~FreeFormInversionCommunicator() { }

protected:	
	errut::bool_t runGA_next(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 const std::shared_ptr<eatk::FitnessComparison> &comparison,
						 int popSize, const std::string &lensFitnessObjectType,
						 const std::string &calculatorType,
	                     grale::LensGACalculatorFactory &calcFactory, 
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator0,
						 const std::vector<uint8_t> &factoryParamBytes,
						 const std::vector<std::unique_ptr<grale::EAParameters>> &allEAParams,
						 const std::vector<grale::LensGAConvergenceParameters> &allConvParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::vector<std::string> &allEATypes) override
	{
		// TODO: Check type name and parameters
		for (size_t i = 0 ; i < allEATypes.size() ; i++)
		{
			// TODO
		}

		std::shared_ptr<grale::LensInversionGAFactoryCommon> genomeCalculator = std::dynamic_pointer_cast<grale::LensInversionGAFactoryCommon>(genomeCalculator0);
		if (!genomeCalculator.get())
			return "Calculator does not seem to be of a type derived from LensInversionGAFactoryCommon";

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

			creation = std::move(lensGACreation);

			if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
									factoryParamBytes, *creation, calc)))
				return "Can't get calculator: " + r.getErrorString();
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
			else if (eaType == "NSGA2" || eaType == "NSGA2-X")
				std::tie(r, previousBestIndividuals, reuseCreation, numGen) = runGA_NSGA2(rng, *creation, comparison, *calc, popSize, genomeCalculator->allowNegativeValues(), genomeCalculator->getNumberOfObjectives(),
								*(allEAParams[i]), allConvParams[i], multiPopParams, eaType, generationCount, previousBestIndividuals);
			else
				r = "Unknown EA type '" + eaType + "', should be either GA, DE or JADE";

			if (!r)
				break;

			std::swap(creation, reuseCreation); // Make sure the next algorithm starts where this one left off

			generationCount += numGen;
		}

		// Note: m_best must be set inside the subroutines; previousBest can be the one from multiple
		//       populations, don't want to recalculate the non-dominated set here

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
						 const std::shared_ptr<eatk::FitnessComparison> &comparison,
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
			WriteLineStdout("GAMESSAGESTR:WARNING: not using previous best for initialization of DE/JADE");

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
						 const std::shared_ptr<eatk::FitnessComparison> &comparison,
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

	std::tuple<errut::bool_t,
		       std::vector<std::vector<std::shared_ptr<eatk::Individual>>>,
			   std::unique_ptr<eatk::IndividualCreation>,
			   size_t> runGA_NSGA2(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 eatk::IndividualCreation &creation,
						 const std::shared_ptr<eatk::FitnessComparison> &comparison,
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
			WriteLineStdout("GAMESSAGESTR:WARNING: not using previous best for initialization of NSGA2");
		if (multiPopParams.get())
			return E("NSGA2 only works with a single population");

		errut::bool_t r;

		Stop stop(-1, generationOffsetForReporting);
		if (!(r = stop.initialize(numObjectives, convParams)))
			return E("Can't initialize stop criterion: " + r.getErrorString());

		std::shared_ptr<eatk::GenomeCrossover> crossOver;
		std::shared_ptr<eatk::GenomeMutation> mutation;

		if (dynamic_cast<const grale::NSGA2Parameters*>(&eaParams))
		{
			const grale::NSGA2Parameters &params = static_cast<const grale::NSGA2Parameters&>(eaParams);
			float mutAmp = params.getSmallMutationSize();
			bool absMut = (mutAmp > 0)?false:true;

			crossOver = std::make_shared<grale::LensGAGenomeCrossover>(rng, allowNegative);
			mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 1.0f, allowNegative, mutAmp, absMut);

			WriteLineStdout("GAMESSAGESTR:Running NSGA2 algorithm, small mutation size = " + std::to_string(params.getSmallMutationSize()));
		}
		else if (dynamic_cast<const grale::NSGA2DELikeCrossoverParameters*>(&eaParams))
		{
			const grale::NSGA2DELikeCrossoverParameters &params = static_cast<const grale::NSGA2DELikeCrossoverParameters&>(eaParams);
			bool extraParent = params.useExtraParent();
			float F = params.getF();
			float CR = params.getCR();
			if (allowNegative)
			{
				crossOver = std::make_shared<grale::LensGAGenomeDELikeCrossover>(rng, extraParent, F, CR);
			}
			else
			{
				// Set lower bound, we need to know the number of parameters in
				// a genome
				auto g = creation.createUnInitializedGenome();
				std::vector<float> lowerBound(grale::LensGAGenome::getSize(*g), 0.0f);

				crossOver = std::make_shared<grale::LensGAGenomeDELikeCrossover>(rng, extraParent, F, CR, lowerBound);
			}

			// No separate mutation, crossover operator does all DE-like operations
			mutation = nullptr;

			std::string Fstr = (std::isnan(F))?std::string("random"):std::to_string(F);
			std::string CRstr = (std::isnan(CR))?std::string("random"):std::to_string(CR);
			WriteLineStdout("GAMESSAGESTR:Running NSGA2 algorithm with experimental DE-like crossover, useextraparent = "
			                + std::to_string(extraParent) + ", F = " + Fstr + ", CR = " + CRstr);
		}
		else
			return E("Invalid EA parameters for NSGA2");

		std::unique_ptr<eatk::PopulationEvolver> evolver = std::make_unique<grale::LensNSGA2Evolver>(rng, crossOver, mutation, comparison, numObjectives);

		MyGA ga;
		if (!(r = ga.run(creation, *evolver, calc, stop, popSize, popSize, popSize*2)))
			return E("Error running NSGA2: " + r.getErrorString());

		std::vector<std::vector<std::shared_ptr<eatk::Individual>>> allBest = { { evolver->getBestIndividuals() } };
		m_best = evolver->getBestIndividuals();

		return { true, allBest,
				 std::make_unique<ReuseCreation>(ga.getPopulations()),
				 ga.getNumberOfGenerations() };
	}
};
