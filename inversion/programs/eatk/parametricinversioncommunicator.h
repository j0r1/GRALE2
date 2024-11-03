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

	~MyGW()
	{
		std::stringstream ss;
		ss << "GAMESSAGESTR: wrote " << m_generationsWritten << " generations of " << m_numWalkers << " walkers, total of " << m_samplesWritten << " samples, each " << m_floatsPerSample << " float values";
		WriteLineStdout(ss.str());
	}

	errut::bool_t init(const std::string &samplesFileName, size_t burninGenerations, size_t annealScale,
			    double alpha0, double alphaMax)
	{
		m_filename = samplesFileName;
		if (!m_filename.empty())
		{
			m_sampleFile.open(m_filename, std::ios::out | std::ios::binary);
			if (!m_sampleFile.is_open())
				return "Unable to open file '" + m_filename + "'";
		}

		m_burnInGenerations = burninGenerations;
		m_annealScaleGenerations = annealScale;
		m_alpha0 = alpha0;
		m_alphaMax = alphaMax;
		m_currentAlpha = 1.0;

		m_generationCount = 0;
		updateAnnealFactor();
		return true;
	}

	void updateAnnealFactor()
	{
		if (m_annealScaleGenerations == 0)
			return;

		double f = (double)m_generationCount/(double)m_annealScaleGenerations;
		double alpha = m_alpha0*pow(1.0/m_alpha0, f);
		if (alpha > m_alphaMax)
			alpha = m_alphaMax;

		m_currentAlpha = alpha;
		setAnnealingExponent(alpha);
	}

	void onSamples(const std::vector<std::shared_ptr<eatk::Individual>> &samples)
	{           
		std::stringstream ss;
		ss << "Generation " << m_generationCount << ": best = ";

		const auto &best = getBestIndividuals();
		if (best.size() == 0)
			ss << "(no best yet)";
		else
			ss << best[0]->fitness()->toString();

		ss << " first = " << samples[0]->fitness()->toString();
		double worst = samples[0]->fitness()->getRealValue(0); // TODO: other objective?
		size_t worstIdx = 0;
		for (size_t i = 0 ; i < samples.size() ; i++)
		{
			double v = samples[i]->fitness()->getRealValue(0); // TODO: other objective?
			if (v > worst)
			{
				worst = v;
				worstIdx = i;
			}
		}
		ss << " worst = " << samples[worstIdx]->fitness()->toString();

		if (m_generationCount < m_burnInGenerations)
			ss << " (burn-in)";
		if (m_currentAlpha != 1.0)
			ss << " (annealing exponent: " << m_currentAlpha << ")";

		WriteLineStdout("GAMESSAGESTR:" + ss.str());

		if (m_generationCount >= m_burnInGenerations && m_sampleFile.is_open())
		{
			m_numWalkers = samples.size();
			for (auto &i : samples)
			{
				assert(dynamic_cast<const eatk::FloatVectorGenome *>(i->genomePtr()));
				const eatk::FloatVectorGenome &genome = static_cast<const eatk::FloatVectorGenome&>(i->genomeRef());
				const std::vector<float> &values = genome.getValues();

				m_sampleFile.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(float));
				m_floatsPerSample = values.size();
				m_samplesWritten++;
			}
			m_generationsWritten++;
		}

		m_generationCount++;
		updateAnnealFactor();
	}
private:
	size_t m_generationCount = 0;
	std::string m_filename;
	std::ofstream m_sampleFile;
	size_t m_burnInGenerations;
	size_t m_annealScaleGenerations;
	double m_alpha0;
	double m_alphaMax;
	double m_currentAlpha;
	size_t m_floatsPerSample = 0;
	size_t m_samplesWritten = 0;
	size_t m_generationsWritten = 0;
	size_t m_numWalkers = 0;
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
		int annealScale = mcmcParams.getAnnealGenerationsTimeScale(); // These don't need *2, is used when a full sample is received
		int burninGenerations = mcmcParams.getBurnInGenerations();
		double a = mcmcParams.getGoodmanWeare_a();
		double alpha0 = mcmcParams.getAnnealAlpha0();
		double alphaMax = mcmcParams.getAnnealAlphaMax();
		std::string samplesFileName = mcmcParams.getSamplesFilename();

		if (maxEAGenerations == 0)
			return "No number of generations to sample has been set";
		if (samplesFileName.empty())
		{
			std::cerr << "\n\n\nWARNING: will not be writing any samples to a file!\n\n\n\n";
			std::this_thread::sleep_for(std::chrono::seconds(3));
		}

		eatk::FixedGenerationsStopCriterion stop(maxEAGenerations);
		MyGW gw(rng, eatk::GoodmanWeareEvolver::NegativeLog, a);
		errut::bool_t r;
		if (!(r = gw.init(samplesFileName, burninGenerations, annealScale, alpha0, alphaMax)))
			return r;

		MyGA ea;
		r = ea.run(creation, gw, *popCalc, stop, popSize, popSize, popSize*3/2);
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
