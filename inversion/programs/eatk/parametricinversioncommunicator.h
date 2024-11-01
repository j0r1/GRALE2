#pragma once

#include "inversioncommunicatorbase.h"
#include "lensgaparametricsingleplanecalculator.h"
#include "mcmcparameters.h"
#include <eatk/vectorgenomefitness.h>
#include <eatk/vectordifferentialevolution.h>
#include <eatk/goodmanweareevolver.h>

class ParametricCreation : public eatk::VectorDifferentialEvolutionIndividualCreation<float,float>
{
public:
	ParametricCreation(size_t numObjectives,
					   const std::vector<float> &lower, const std::vector<float> &upper,
	                   const std::shared_ptr<eatk::RandomNumberGenerator> &rng)
		: eatk::VectorDifferentialEvolutionIndividualCreation<float,float>(lower, upper, rng),
		  m_numObj(numObjectives)
	{
	}
	
	std::shared_ptr<eatk::Fitness> createEmptyFitness() override
	{
		return std::make_shared<eatk::FloatVectorFitness>(m_numObj);
	}

	size_t m_numObj;
};

class MyGW : public eatk::GoodmanWeareEvolver
{       
public: 
	MyGW(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, ProbType probType, double a)
		: eatk::GoodmanWeareEvolver(rng, probType, a) { } 

	void onSamples(const std::vector<std::shared_ptr<eatk::Individual>> &samples)
	{           
		std::stringstream ss;
		ss << "Generation " << generationCount++ << ": best = ";

		const auto &best = getBestIndividuals();
		if (best.size() == 0)
			ss << "(no best yet)";
		else
			ss << best[0]->fitness()->toString();
		ss << " first walker = " << samples[0]->fitness()->toString();
		WriteLineStdout("GAMESSAGESTR:" + ss.str());
		// TODO: log the samples!
	}
private:
	size_t generationCount = 0;
};          

class ParametricInversionCommunicator : public InversionCommunicatorBase
{
public:
	ParametricInversionCommunicator() { }
	~ParametricInversionCommunicator() { }

	static std::shared_ptr<eatk::Genome> getReferenceGenomeForMPI()
	{
		return std::make_shared<eatk::FloatVectorGenome>();
	}

	static std::shared_ptr<eatk::Fitness> getReferenceFitnessForMPI()
	{
		return std::make_shared<eatk::FloatVectorFitness>();
	}


protected:	

	errut::bool_t runMCMC(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						int popSize, const std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> &calculator,
						const std::shared_ptr<eatk::PopulationFitnessCalculation> &popCalc,
						eatk::IndividualCreation &creation,
						const grale::EAParameters &eaParams)

	{
		if (calculator->getNumberOfObjectives() != 1)
			return "Only a single objective can be used for MCMC";

		const grale::LensFitnessObject &fitObj = calculator->getFitnessObject();
		if (!fitObj.isNegativeLogProb_Overall(0))
			return "The fitness measure is not the negative of a (log) probability";

		if (!dynamic_cast<const grale::MCMCParameters*>(&eaParams))
			return "Parameters are not of type MCMCParameters";
		const grale::MCMCParameters &mcmcParams = static_cast<const grale::MCMCParameters&>(eaParams);

		// NOTE: we're ignoring the convergence parameters!
		
		int maxEAGenerations = mcmcParams.getSampleGenerations() * 2; // Here we use *2 since it takes 2 generations to advance the full set of walkers
		if (maxEAGenerations == 0)
			return "No number of generations to sample has been set";

		eatk::FixedGenerationsStopCriterion stop(maxEAGenerations);

		double a = mcmcParams.getGoodmanWeare_a();
		std::string samplesFileName = mcmcParams.getSamplesFilename();
		if (samplesFileName.empty())
		{
			std::cerr << "\n\n\nWARNING: will not be writing any samples to a file!\n\n\n\n";
			std::this_thread::sleep_for(std::chrono::seconds(3));
		}
		MyGW gw(rng, eatk::GoodmanWeareEvolver::NegativeLog, a); // TODO: make a file to log the samples

		MyGA ea;
		bool_t r = ea.run(creation, gw, *popCalc, stop, popSize, popSize, popSize*3/2);
		if (!r)
			return r;

		m_best = gw.getBestIndividuals();
		return true;
	}

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

		if (allEATypes.size() != 1)
			return "Currently only a single EA type can be used for parametric inversion";

		if (multiPopParams.get())
			return "Parametric inversion only works with a single population";

		std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> genomeCalculator = std::dynamic_pointer_cast<grale::LensGAParametricSinglePlaneCalculator>(genomeCalculator0);
		if (!genomeCalculator)
			return "Calculator does not seem to be a LensGAParametricSinglePlaneCalculator one";

		size_t numObjectives = genomeCalculator->getNumberOfObjectives();

		std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
		errut::bool_t r;

		std::unique_ptr<eatk::IndividualCreation> creation = std::make_unique<ParametricCreation>(numObjectives, genomeCalculator->getInitMin(), genomeCalculator->getInitMax(), rng);

		if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
									factoryParamBytes, *creation, calc)))
			return "Can't get calculator: " + r.getErrorString();

		const grale::EAParameters &eaParams = *allEAParams[0];
		const grale::LensGAConvergenceParameters &convParams = allConvParams[0];
		std::string eaType = allEATypes[0];

		// Do MCMC in another routine
		if (eaType == "MCMC")
			return runMCMC(rng, popSize, genomeCalculator, calc, *creation, eaParams);

		// Continue with JADE or DE
		if (!(eaType == "JADE" || eaType == "DE"))
			return "Currently only MCMC, DE or JADE is supported";

		Stop stop(-1, 0);
		if (!(r = stop.initialize(numObjectives, convParams)))
			return "Can't initialize stop criterion: " + r.getErrorString();
		
		auto mut = std::make_shared<eatk::VectorDifferentialEvolutionMutation<float>>();
		auto cross = std::make_shared<eatk::VectorDifferentialEvolutionCrossover<float>>(rng, genomeCalculator->getHardMin(), genomeCalculator->getHardMax());

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
				evolver = std::make_unique<eatk::JADEEvolver>(rng, mut, cross, comparison, 0,
						                                           p, c, useArch, initMuF, initMuCR,
																   1, nullptr, needStrictlyBetter);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<eatk::JADEEvolver>(rng, mut, cross, comparison,
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
				evolver = std::make_unique<eatk::DifferentialEvolutionEvolver>(rng, mut, F, cross, CR, comparison,
						                                         0, 1, nullptr, needStrictlyBetter);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<eatk::DifferentialEvolutionEvolver>(rng, mut, params.getF(), cross, params.getCR(), comparison,
						                                         -1, numObjectives, ndCreator,
																 needStrictlyBetter); // -1 signals multi-objective
			}
		}
		else
			return "Unexpected eaType '" + eaType + "'";

		MyGA ga;
		if (!(r = ga.run(*creation, *evolver, *calc, stop, popSize, popSize, popSize*2)))
			return "Error running GA: " + r.getErrorString();

		m_best = evolver->getBestIndividuals();

		return true;
	}
};
