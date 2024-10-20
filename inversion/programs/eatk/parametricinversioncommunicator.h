#pragma once

#include "inversioncommunicatorbase.h"
#include "lensgaparametricsingleplanecalculator.h"

class ParametricInversionCommunicator : public InversionCommunicatorBase
{
public:
	ParametricInversionCommunicator() { }
	~ParametricInversionCommunicator() { }

protected:	

	errut::bool_t runGA(int popSize, const std::string &lensFitnessObjectType,
						 const std::string &calculatorType,
	                     grale::LensGACalculatorFactory &calcFactory, 
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator0,
						 const std::vector<uint8_t> &factoryParamBytes,
						 const std::vector<std::unique_ptr<grale::EAParameters>> &allEAParams,
						 const std::vector<grale::LensGAConvergenceParameters> &allConvParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::vector<std::string> &allEATypes)
	{
		return "TODO: parametric inversion is under construction";

		if (allEATypes.size() != allEAParams.size() || allEATypes.size() != allConvParams.size())
			return "Unexpected mismatch between number of EA types (" + std::to_string(allEATypes.size()) +
				   "), EA parameters (" + std::to_string(allEAParams.size()) + ") and convergence parameters (" + 
				   std::to_string(allConvParams.size()) + ")";

		// TODO: Check type name and parameters
		for (size_t i = 0 ; i < allEATypes.size() ; i++)
		{
			// TODO
		}

		std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> genomeCalculator = std::dynamic_pointer_cast<grale::LensGAParametricSinglePlaneCalculator>(genomeCalculator0);
		if (!genomeCalculator.get())
			return "Calculator does not seem to be a LensGAParametricSinglePlaneCalculator one";

		std::shared_ptr<grale::RandomNumberGenerator> rng0 = std::make_shared<grale::RandomNumberGenerator>();
		WriteLineStdout("GAMESSAGESTR:RNG SEED: " + std::to_string(rng0->getSeed()));

#if 0
		std::shared_ptr<eatk::RandomNumberGenerator> rng = std::make_shared<RngWrapper>(rng0);
#else
		std::shared_ptr<eatk::RandomNumberGenerator> rng = rng0;
#endif

		/*
		auto comparison = std::make_shared<grale::LensGAFitnessComparison>();
		m_selector = std::make_shared<SubsequentBestIndividualSelector>(
								genomeCalculator->getNumberOfObjectives(),
								comparison);
		*/
		// After this is created, calculatorCleanup() should be called as well (for MPI at the moment)
		std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
		errut::bool_t r;

		/*
		std::unique_ptr<eatk::IndividualCreation> creation;
		{
			std::unique_ptr<grale::LensGAIndividualCreation> lensGACreation = std::make_unique<grale::LensGAIndividualCreation>(rng, 
							  genomeCalculator->getNumberOfBasisFunctions(),
							  genomeCalculator->getNumberOfSheets(),
							  genomeCalculator->allowNegativeValues(),
							  genomeCalculator->getNumberOfObjectives());

			if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
									factoryParamBytes, *creation, calc)))
				return "Can't get calculator: " + r.getErrorString();

			creation = std::move(lensGACreation);
		}
		*/

		// TODO: Call DE or JADE code

		// Note: m_best must be set inside the subroutines; previousBest can be the one from multiple
		//       populations, don't want to recalculate the non-dominated set here

		calculatorCleanup();
		return r;
	}

	/*

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
*/
};
