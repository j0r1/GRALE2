#pragma once

#include "inversioncommunicatorbase.h"
#include "lensgaparametricsingleplanecalculator.h"
#include "mcmcparameters.h"
#include <eatk/vectorgenomefitness.h>
#include <eatk/vectordifferentialevolution.h>
#include <eatk/goodmanweareevolver.h>
#include <eatk/metropolishastingsevolver.h>
#include <eatk/vectorgenomeuniformmutation.h>
#include <eatk/vectorgenomefractionalmutation.h>
#include <eatk/vectorgenomeuniformcrossover.h>
#include <eatk/vectorgenomedelikecrossover.h>
#include <eatk/populationreusecreation.h>
#include <eatk/singlepopulationcrossover.h>
#include <eatk/rankparentselection.h>
#include <eatk/simplesortedpopulation.h>
#include <eatk/singlebestelitism.h>
#include <eatk/remainingtargetpopulationsizeiteration.h>
#include <iomanip>
#include <list>

class DummyEvolver : public eatk::PopulationEvolver
{
public:
	DummyEvolver(size_t expectedPopSize) : m_expectedPopSize(expectedPopSize) { }

	errut::bool_t check(const std::shared_ptr<eatk::Population> &population)
	{
		if (population->size() != m_expectedPopSize)
			return "Expecting a population size of " + std::to_string(m_expectedPopSize);
		return true;
	}

	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize)
	{
		if (generation != 0)
			return "Expecting only a single generation";
		if (population->size() != m_expectedPopSize)
			return "Expecting a population size of " + std::to_string(m_expectedPopSize);
		if (targetPopulationSize != m_expectedPopSize)
			return "Expecting a target population size of " + std::to_string(m_expectedPopSize);

		m_best.clear();

		// We'll store all members in the best individuals. This is actually
		// meant for the non-dominated set, but we can abuse this to return the
		// fitness values for a whole set of genomes
		for (auto &i : population->individuals())
		{
			m_best.push_back(i->createCopy());

			// For now, an assertion is triggered if the OpenCL code doesn't
			// need to perform any calculations. This currently makes sure
			// that a recalculation is needed
			i->fitness()->setCalculated(false);
		}

		return true;
	}

	const std::vector<std::shared_ptr<eatk::Individual>> &getBestIndividuals() const { return m_best; }
private:
	std::vector<std::shared_ptr<eatk::Individual>> m_best;
	size_t m_expectedPopSize;
};

class CustomFloatVectorFitness : public eatk::FloatVectorFitness
{
public:
	CustomFloatVectorFitness(size_t n = 0) : eatk::FloatVectorFitness(n) { }

	// Just so that higher precision is used when reporting fitness.
	// When it is parsed from the string, the difference with the actual
	// value will be smaller.
	std::string toString() const override
	{
		if (!Fitness::isCalculated())
			return "?";

		std::stringstream ss;

		ss << "[";
		for (auto x : m_values)
			ss << " " << std::setprecision(15) << x;
		ss << " ]";

		return ss.str();
	}

	std::shared_ptr<Fitness> createCopy(bool copyContents = true) const override
	{
		auto g = ValueVector<Fitness, float>::template createCopy<CustomFloatVectorFitness>(copyContents);
		if (copyContents && Fitness::isCalculated())
			g->setCalculated();
		return g;
	}
};

class PredefinedIndividualCreation : public eatk::IndividualCreation
{
public:
	PredefinedIndividualCreation(const std::vector<std::vector<float>> &genomes, size_t genomeSize, size_t fitnessSize)
	 : m_genomeSize(genomeSize), m_fitnessSize(fitnessSize)
	{
		for (auto &g : genomes)
			m_genomes.push_back(g);
	}

	std::shared_ptr<eatk::Genome> createInitializedGenome() override
	{
		if (m_genomes.empty())
			return nullptr;

		std::vector<float> values = m_genomes.front();
		m_genomes.pop_front();

		auto g = std::make_shared<eatk::FloatVectorGenome>(m_genomeSize);
		std::vector<float> &genomeValues = g->getValues();

		if (values.size() != genomeValues.size())
			return nullptr;

		for (size_t i = 0 ; i < values.size() ; i++)
			genomeValues[i] = values[i];

		return g;
	}

	std::shared_ptr<eatk::Genome> createUnInitializedGenome() override
	{
		return std::make_shared<eatk::FloatVectorGenome>(m_genomeSize);
	}

	std::shared_ptr<eatk::Fitness> createEmptyFitness() override
	{
		return std::make_shared<CustomFloatVectorFitness>(m_fitnessSize);
	}
private:
	size_t m_genomeSize, m_fitnessSize;
	std::list<std::vector<float>> m_genomes;
};

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
		return std::make_shared<CustomFloatVectorFitness>(m_numObj);
	}

	size_t m_numObj;
};

template <class BaseClass>
class CommonSamplingCode : public BaseClass
{
public:
	template<typename... Args>
	CommonSamplingCode(Args&&... args) : BaseClass(std::forward<Args>(args)...) { }

	~CommonSamplingCode()
	{
		std::stringstream ss;
		ss << "GAMESSAGESTR: wrote " << m_generationsWritten << " generations of " << m_numWalkers << " walkers, total of " << m_samplesWritten << " samples, each " << m_floatsPerSample << " float values";
		WriteLineStdout(ss.str());
	}

	errut::bool_t init(const std::string &samplesFileName, 
			    const std::string &logProbFn,
			    size_t burninGenerations, size_t annealScale,
			    double alpha0, double alphaMax)
	{
		m_filename = samplesFileName;
		m_logProbFn = logProbFn;
		if (!m_filename.empty())
		{
			m_sampleFile.open(m_filename, std::ios::out | std::ios::binary);
			if (!m_sampleFile.is_open())
				return "Unable to open samples file '" + m_filename + "'";
		}
		if (!m_logProbFn.empty())
		{
			m_logProbFile.open(m_logProbFn, std::ios::out | std::ios::binary);
			if (!m_logProbFile.is_open())
				return "Unable to open log prob file '" + m_logProbFn + "'";
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
		BaseClass::setAnnealingExponent(alpha);
	}

	void onSamples(const std::vector<std::shared_ptr<eatk::Individual>> &samples)
	{           
		std::stringstream ss;
		ss << "Generation " << m_generationCount << ": best = ";

		// Note: when sampling it is possible that this best performing genome
		//       is not part of recorded samples. It may be that the best one
		//       was encountered during the burn-in phase. Also, if the
		//       Goodman-Weare sampler is used the best genome is not necessarily
		//       accepted into the chain. It's still saved as the best though.

		const auto &best = BaseClass::getBestIndividuals();
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

		if (m_generationCount >= m_burnInGenerations && (m_sampleFile.is_open() || m_logProbFile.is_open()))
		{
			// TODO If origin parameters are used, transform these to actual parameters?
			m_numWalkers = samples.size();

			if (m_sampleFile.is_open())
			{
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

			if (m_logProbFile.is_open())
			{
				m_logProbBuffer.clear();

				for (auto &i : samples)
				{
					double v = i->fitness()->getRealValue(0); // TODO: other objective?
					m_logProbBuffer.push_back(-(float)v); // Write actual log probs, fitness is the negative
				}

				m_logProbFile.write(reinterpret_cast<const char*>(m_logProbBuffer.data()), m_logProbBuffer.size() * sizeof(float));
			}
		}

		m_generationCount++;
		updateAnnealFactor();

		// if (best.size() > 0)
		// 	std::cerr << "    BEST: " << best[0]->toString() << std::endl;
		// for (auto &s : samples)
		// 	std::cerr << "   " << s->toString() << std::endl;
	}
private:
	size_t m_generationCount = 0;
	std::string m_filename, m_logProbFn;
	std::ofstream m_sampleFile, m_logProbFile;
	size_t m_burnInGenerations = 0;
	size_t m_annealScaleGenerations = 0;
	double m_alpha0 = 1.0;
	double m_alphaMax = 1.0;
	double m_currentAlpha = 1.0;
	size_t m_floatsPerSample = 0;
	size_t m_samplesWritten = 0;
	size_t m_generationsWritten = 0;
	size_t m_numWalkers = 0;

	std::vector<float> m_logProbBuffer;
};

class MyGW : public CommonSamplingCode<eatk::GoodmanWeareEvolver>
{       
public: 
	MyGW(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, ProbType probType, double a)
		: CommonSamplingCode<eatk::GoodmanWeareEvolver>(rng, probType, a) { } 
};          

class MyMH : public CommonSamplingCode<eatk::MetropolisHastingsEvolver>
{
public:
	MyMH(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
		 const std::vector<double> &stepScales,
		 ProbType probType)
		: CommonSamplingCode<eatk::MetropolisHastingsEvolver>(rng, stepScales, probType) { }
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
		return std::make_shared<CustomFloatVectorFitness>();
	}


protected:	

	errut::bool_t runMCMC(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						int popSize, const std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> &calculator,
						const std::shared_ptr<eatk::PopulationFitnessCalculation> &popCalc,
						std::shared_ptr<eatk::IndividualCreation> &creation, const std::string &eaType,
						const grale::EAParameters &eaParams)

	{
		if (calculator->getNumberOfObjectives() != 1)
			return "Only a single objective can be used for MCMC";

		const grale::LensFitnessObject &fitObj = calculator->getFitnessObject();
		if (!fitObj.isNegativeLogProb_Overall(0))
			return "The fitness measure is not the negative of a (log) probability";

		if (calculator->isRandomizingInputPositions())
			return "MCMC sampling cannot be used together with input image randomization";

		if (!dynamic_cast<const grale::GeneralMCMCParameters*>(&eaParams))
			return "Parameters are not of type GeneralMCMCParameters";
		const grale::GeneralMCMCParameters &mcmcParams = static_cast<const grale::GeneralMCMCParameters&>(eaParams);

		// NOTE: we're ignoring the convergence parameters!
		
		int maxEAGenerations = mcmcParams.getSampleGenerations();
		int annealScale = mcmcParams.getAnnealGenerationsTimeScale(); // These don't need *2, is used when a full sample is received
		int burninGenerations = mcmcParams.getBurnInGenerations();
		double alpha0 = mcmcParams.getAnnealAlpha0();
		double alphaMax = mcmcParams.getAnnealAlphaMax();
		std::string samplesFileName = mcmcParams.getSamplesFilename();
		std::string logProbFn = mcmcParams.getLogProbFilename();

		if (maxEAGenerations == 0)
			return "No number of generations to sample has been set";
		if (samplesFileName.empty())
		{
			std::cerr << "\n\n\nWARNING: will not be writing any samples to a file!\n\n\n\n";
			std::this_thread::sleep_for(std::chrono::seconds(3));
		}

		std::shared_ptr<eatk::PopulationEvolver> evolver;
		errut::bool_t r;
		size_t popMax = popSize;

		if (eaType == "MCMC")
		{
			if (!dynamic_cast<const grale::MCMCParameters*>(&eaParams))
				return "Parameters are not of type MCMCParameters";
			const grale::MCMCParameters &mcmcParams = static_cast<const grale::MCMCParameters&>(eaParams);
			double a = mcmcParams.getGoodmanWeare_a();	

			std::shared_ptr<MyGW> gw = std::make_shared<MyGW>(rng, eatk::SamplingEvolver::NegativeLog, a);
			if (!(r = gw->init(samplesFileName, logProbFn, burninGenerations, annealScale, alpha0, alphaMax)))
				return r;

			evolver = gw;
			maxEAGenerations *= 2; // Here we use *2 since it takes 2 generations to advance the full set of walkers
			popMax = popSize*3/2;
			WriteLineStdout("GAMESSAGESTR: Running MCMC algorithm"); // TODO: more info about settings
		}
		else if (eaType == "MCMC-MH")
		{
			if (!dynamic_cast<const grale::MetropolisHastingsMCMCParameters*>(&eaParams))
				return "Parameters are not of type MetropolisHastingsMCMCParameters";
			const grale::MetropolisHastingsMCMCParameters &mcmcParams = static_cast<const grale::MetropolisHastingsMCMCParameters&>(eaParams);
			const std::vector<double> stepScales = mcmcParams.getStepScales();

			std::shared_ptr<MyMH> mh = std::make_shared<MyMH>(rng, stepScales, eatk::SamplingEvolver::NegativeLog);
			if (!(r = mh->init(samplesFileName, logProbFn, burninGenerations, annealScale, alpha0, alphaMax)))
				return r;

			evolver = mh;
			popMax = popSize*2;
			WriteLineStdout("GAMESSAGESTR: Running MCMC-MH algorithm"); // TODO: more info about settings
		}
		else
			return "Unrecognized MCMC EA type '" + eaType + "'";

		assert(evolver.get());

		eatk::FixedGenerationsStopCriterion stop(maxEAGenerations);
		MyGA ea;
		r = ea.run(*creation, *evolver, *popCalc, stop, popSize, popSize, popMax);
		if (!r)
			return r;

		m_best = evolver->getBestIndividuals();
		creation = std::make_shared<eatk::PopulationReuseCreation>(ea.getPopulations());
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
		if (multiPopParams.get())
			return "Parametric inversion only works with a single population";

		std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> genomeCalculator = std::dynamic_pointer_cast<grale::LensGAParametricSinglePlaneCalculator>(genomeCalculator0);
		if (!genomeCalculator)
			return "Calculator does not seem to be a LensGAParametricSinglePlaneCalculator one";

		size_t numObjectives = genomeCalculator->getNumberOfObjectives();

		std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
		errut::bool_t r;

		std::shared_ptr<eatk::IndividualCreation> creation = std::make_unique<ParametricCreation>(numObjectives, genomeCalculator->getInitMin(), genomeCalculator->getInitMax(), rng);

		if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
									factoryParamBytes, *creation, calc)))
			return "Can't get calculator: " + r.getErrorString();

		if (allEATypes.size() == 1 && allEATypes[0] == "CALCULATE")
		{
			if (!(r = runCalculate(popSize, genomeCalculator, calc, creation)))
				return r;
			return true;
		}

		if (genomeCalculator->getGenomesToCalculateFitness().size() > 0)
			return "Genomes to calculate fitness for is not allowed to be used for EA types other than 'CALCULATE'";
		// TODO: Check type name and parameters beforehand?

		for (size_t i = 0 ; i < allEATypes.size() ; i++)
		{
			const grale::EAParameters &eaParams = *allEAParams[i];
			const grale::LensGAConvergenceParameters &convParams = allConvParams[i];
			std::string eaType = allEATypes[i];
			bool_t r;

			// Do MCMC in another routine
			if (eaType == "MCMC" || eaType == "MCMC-MH")
			{
				if (!(r = runMCMC(rng, popSize, genomeCalculator, calc, creation, eaType, eaParams)))
					return r;
			}
			else
			{
				Stop stop(-1, 0);
				if (!(r = stop.initialize(numObjectives, convParams)))
					return "Can't initialize stop criterion: " + r.getErrorString();

				if (eaType == "JADE" || eaType == "DE")
				{
					if (!(r = runDE(eaType, numObjectives, rng, comparison, popSize, genomeCalculator, calc, creation, eaParams, stop)))
						return r;
				}
				else if (eaType == "NSGA2" || eaType == "NSGA2-X")
				{
					if (!(r = runNSGA2(eaType, numObjectives, rng, comparison, popSize, genomeCalculator, calc, creation, eaParams, stop)))
						return r;
				}
				else if (eaType == "GA")
				{
					if (!(r = runGA(eaType, numObjectives, rng, comparison, popSize, genomeCalculator, calc, creation, eaParams, stop)))
						return r;
				}
				else
					return "Unexpected EA type '" + eaType + "'";
			}
		}

		return true;
	}

	errut::bool_t runCalculate(int popSize, const std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> &genomeCalculator,
						const std::shared_ptr<eatk::PopulationFitnessCalculation> &calc,
						std::shared_ptr<eatk::IndividualCreation> &creation)
	{
		eatk::FixedGenerationsStopCriterion stop(0);
		MyGA ea;
		size_t numGenomesToCalculate = genomeCalculator->getGenomesToCalculateFitness().size();

		if (numGenomesToCalculate == 0)
		{
			// In this case, the lens model limits are set in such a way
			// that only one specific lens can be generated

			DummyEvolver evolver(1); // Just to calculate fitness values
			bool_t r = ea.run(*creation, evolver, *calc, stop, popSize, popSize, popSize);
			if (!r)
				return r;

			m_best = evolver.getBestIndividuals();
		}
		else
		{
			std::shared_ptr<eatk::Genome> refGenome = creation->createUnInitializedGenome();
			std::shared_ptr<eatk::Fitness> refFitness = creation->createEmptyFitness();
			eatk::FloatVectorGenome *pRefFloatGenome = dynamic_cast<eatk::FloatVectorGenome*>(refGenome.get());
			eatk::FloatVectorFitness *pRefFloatFitness = dynamic_cast<eatk::FloatVectorFitness*>(refFitness.get());
			if (pRefFloatGenome == nullptr)
				return "Error in fitness calculation check: reference genome is not of type FloatVectorGenome";
			if (pRefFloatFitness == nullptr)
				return "Error in fitness calculation check: reference fitness is not of type FloatVectorFitness";

			size_t genomeSize = pRefFloatGenome->getValues().size();
			size_t fitnessSize = pRefFloatFitness->getValues().size();

			DummyEvolver evolver(numGenomesToCalculate);
			PredefinedIndividualCreation predefCreation(genomeCalculator->getGenomesToCalculateFitness(),
														genomeSize, fitnessSize);

			bool_t r = ea.run(predefCreation, evolver, *calc, stop, popSize, popSize, popSize);
			if (!r)
				return r;

			m_best = evolver.getBestIndividuals();
		}

		return true;
	}

	errut::bool_t runGA(const std::string &eaType, size_t numObjectives, const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						const std::shared_ptr<eatk::FitnessComparison> &comparison,
						int popSize, const std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> &genomeCalculator,
						const std::shared_ptr<eatk::PopulationFitnessCalculation> &calc,
						std::shared_ptr<eatk::IndividualCreation> &creation,
						const grale::EAParameters &eaParams, eatk::StopCriterion &stop)
	{
		const std::vector<float> &initMin = genomeCalculator->getInitMin();
		const std::vector<float> &initMax = genomeCalculator->getInitMax();
		const std::vector<float> &hardMin = genomeCalculator->getHardMin();
		const std::vector<float> &hardMax = genomeCalculator->getHardMax();

		if (eaType != "GA")
			return "Unexpected EA type '" + eaType + "'";

		if (numObjectives != 1)
			return "Regular GA for parametric inversion only supports a single objective; you can use NSGA2";

		if (!dynamic_cast<const grale::GAParameters*>(&eaParams))
			return "Parameters are not suitable for GA";

		const grale::GAParameters &params = static_cast<const grale::GAParameters&>(eaParams);
		double selPressure = params.getSelectionPressure();
		bool useElitism = params.getUseElitism();
		bool alwaysIncludeBest = params.getAlwaysIncludeBest();
		double crossRate = params.getCrossOverRate();
		float mutAmp = params.getSmallMutationSize();
		bool absMut = (mutAmp > 0)?false:true;

		assert(crossRate >= 0 && crossRate <= 1.0);

		std::shared_ptr<eatk::GenomeMutation> mut;
		// Get the right mutation operator
		{
			auto genome = creation->createUnInitializedGenome();
			assert(dynamic_cast<const eatk::FloatVectorGenome*>(genome.get()));
			const eatk::FloatVectorGenome &vg = static_cast<const eatk::FloatVectorGenome&>(*genome);
			size_t numParameters = vg.getValues().size();
			double mutFrac = 1.0/numParameters;

			if (absMut)
			{
				mut = std::make_shared<eatk::VectorGenomeUniformMutation<float>>(mutFrac, initMin, initMax, rng);
			}
			else
			{
				std::vector<float> refScales(initMin.size());

				// Calculate scale vector
				assert(initMin.size() == initMax.size());
				for (size_t i = 0 ; i < initMin.size() ; i++)
				{
					float diff = std::abs(initMax[i] - initMin[i]);
					refScales[i] = diff*mutAmp;
				}

				struct AdjRoutine
				{
					static float adjustValue(float oldVal, float scale, float hardMin, float hardMax, eatk::RandomNumberGenerator &rng)
					{
						// Same routine as in LensGAGenomeMutation
						float p = rng.getRandomFloat()*2.0f-1.0f;
						float x = std::tan(p*(float)(grale::CONST_PI/4.0))*scale;
						return x + oldVal;
					}
				};

				mut = std::make_shared<eatk::VectorGenomeFractionalMutation<float, AdjRoutine, true>>(mutFrac, refScales, hardMin, hardMax, rng);
			}
		}

		std::shared_ptr<eatk::GenomeCrossover> cross = std::make_shared<eatk::VectorGenomeUniformCrossover<float>>(rng, false);

		std::shared_ptr<eatk::Elitism> elitism;
		if (useElitism)
		{
			if (alwaysIncludeBest)
				elitism = std::make_shared<eatk::SingleBestElitism>(true, mut);
			else
				elitism = std::make_shared<eatk::SingleBestElitism>(false, mut);
		}
		else
		{
			if (alwaysIncludeBest)
				elitism = std::make_shared<eatk::SingleBestElitism>(false, nullptr);
		}

		// TODO: use a wrapper to perform a population check and do the dump population stuff
		std::unique_ptr<eatk::PopulationEvolver> evolver = std::make_unique<eatk::SinglePopulationCrossover>(1.0-crossRate, false,
																											 std::make_shared<eatk::SimpleSortedPopulation>(std::make_shared<eatk::VectorFitnessComparison<float>>()),
																											 std::make_shared<eatk::RankParentSelection>(selPressure, rng),
																											 cross,
																											 mut,
																											 elitism,
																											 std::make_shared<eatk::RemainingTargetPopulationSizeIteration>(),
																											 rng);

		WriteLineStdout("GAMESSAGESTR:Running GA algorithm, selection pressure = " + std::to_string(params.getSelectionPressure()) +
				        ", elitism = " + std::to_string((int)params.getUseElitism()) +
						", always include best = " + std::to_string((int)params.getAlwaysIncludeBest()) + 
						", crossover rate = " + std::to_string(params.getCrossOverRate()) +
						", small mutation size = " + std::to_string(params.getSmallMutationSize()));

		MyGA ga;
		errut::bool_t r;
		if (!(r = ga.run(*creation, *evolver, *calc, stop, popSize, popSize, popSize+2))) // +2 for the elites
			return "Error running GA: " + r.getErrorString();

		m_best = evolver->getBestIndividuals();
		creation = std::make_shared<eatk::PopulationReuseCreation>(ga.getPopulations());

		return true;
	}

	errut::bool_t runNSGA2(const std::string &eaType, size_t numObjectives, const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						const std::shared_ptr<eatk::FitnessComparison> &comparison,
						int popSize, const std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> &genomeCalculator,
						const std::shared_ptr<eatk::PopulationFitnessCalculation> &calc,
						std::shared_ptr<eatk::IndividualCreation> &creation,
						const grale::EAParameters &eaParams, eatk::StopCriterion &stop)
	{
		std::unique_ptr<eatk::PopulationEvolver> evolver;
		std::shared_ptr<eatk::GenomeMutation> mut;
		std::shared_ptr<eatk::GenomeCrossover> cross;

		const std::vector<float> &initMin = genomeCalculator->getInitMin();
		const std::vector<float> &initMax = genomeCalculator->getInitMax();
		const std::vector<float> &hardMin = genomeCalculator->getHardMin();
		const std::vector<float> &hardMax = genomeCalculator->getHardMax();

		if (eaType == "NSGA2")
		{
			if (!dynamic_cast<const grale::NSGA2Parameters*>(&eaParams))
				return "Parameters are not suitable for NSGA2";
		
			const grale::NSGA2Parameters &params = static_cast<const grale::NSGA2Parameters&>(eaParams);
			float mutAmp = params.getSmallMutationSize();
			bool absMut = (mutAmp > 0)?false:true;

			auto genome = creation->createUnInitializedGenome();
			assert(dynamic_cast<const eatk::FloatVectorGenome*>(genome.get()));
			const eatk::FloatVectorGenome &vg = static_cast<const eatk::FloatVectorGenome&>(*genome);
			size_t numParameters = vg.getValues().size();
			double mutFrac = 1.0/numParameters;

			if (absMut)
			{
				mut = std::make_shared<eatk::VectorGenomeUniformMutation<float>>(mutFrac, initMin, initMax, rng);
			}
			else
			{
				std::vector<float> refScales(initMin.size());

				// Calculate scale vector
				assert(initMin.size() == initMax.size());
				for (size_t i = 0 ; i < initMin.size() ; i++)
				{
					float diff = std::abs(initMax[i] - initMin[i]);
					refScales[i] = diff*mutAmp;
				}

				struct AdjRoutine
				{
					static float adjustValue(float oldVal, float scale, float hardMin, float hardMax, eatk::RandomNumberGenerator &rng)
					{
						// Same routine as in LensGAGenomeMutation
						float p = rng.getRandomFloat()*2.0f-1.0f;
						float x = std::tan(p*(float)(grale::CONST_PI/4.0))*scale;
						return x + oldVal;
					}
				};

				mut = std::make_shared<eatk::VectorGenomeFractionalMutation<float, AdjRoutine, true>>(mutFrac, refScales, hardMin, hardMax, rng);
			}

			cross = std::make_shared<eatk::VectorGenomeUniformCrossover<float>>(rng, false);

			WriteLineStdout("GAMESSAGESTR:Running NSGA2 algorithm, small mutation size = " + std::to_string(params.getSmallMutationSize()));
		}
		else if (eaType == "NSGA2-X")
		{
			if (!dynamic_cast<const grale::NSGA2DELikeCrossoverParameters*>(&eaParams))
				return "Parameters are not suitable for NSGA2-X";

			const grale::NSGA2DELikeCrossoverParameters &params = static_cast<const grale::NSGA2DELikeCrossoverParameters&>(eaParams);
			bool extraParent = params.useExtraParent();
			float F = params.getF();
			float CR = params.getCR();

			cross = std::make_shared<eatk::VectorGenomeDELikeCrossOver<float>>(rng, extraParent, F, CR, hardMin, hardMax);

			std::string Fstr = (std::isnan(F))?std::string("random"):std::to_string(F);
			std::string CRstr = (std::isnan(CR))?std::string("random"):std::to_string(CR);
			WriteLineStdout("GAMESSAGESTR:Running NSGA2 algorithm with experimental DE-like crossover, useextraparent = "
			                + std::to_string(extraParent) + ", F = " + Fstr + ", CR = " + CRstr);
		}
		else
			return "Unexpected EA type '" + eaType + "'";

		evolver = std::make_unique<grale::LensNSGA2Evolver>(false, rng, cross, mut, comparison, numObjectives);

		MyGA ga;
		errut::bool_t r;
		if (!(r = ga.run(*creation, *evolver, *calc, stop, popSize, popSize, popSize*2)))
			return "Error running GA: " + r.getErrorString();

		m_best = evolver->getBestIndividuals();
		creation = std::make_shared<eatk::PopulationReuseCreation>(ga.getPopulations());

		return true;
	}

	errut::bool_t runDE(const std::string &eaType, size_t numObjectives, const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						const std::shared_ptr<eatk::FitnessComparison> &comparison,
						int popSize, const std::shared_ptr<grale::LensGAParametricSinglePlaneCalculator> &genomeCalculator,
						const std::shared_ptr<eatk::PopulationFitnessCalculation> &calc,
						std::shared_ptr<eatk::IndividualCreation> &creation,
						const grale::EAParameters &eaParams, eatk::StopCriterion &stop)

	{
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
				evolver = std::make_unique<grale::LensJADEEvolver>(false, rng, mut, cross, comparison, 0,
						                                           p, c, useArch, initMuF, initMuCR,
																   1, nullptr, needStrictlyBetter);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<grale::LensJADEEvolver>(false, rng, mut, cross, comparison,
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
				evolver = std::make_unique<grale::LensDEEvolver>(false, rng, mut, F, cross, CR, comparison,
						                                         0, 1, nullptr, needStrictlyBetter);
			}
			else // multi-objective
			{
				auto ndCreator = std::make_shared<eatk::FasterNonDominatedSetCreator>(comparison, numObjectives);
				evolver = std::make_unique<grale::LensDEEvolver>(false, rng, mut, params.getF(), cross, params.getCR(), comparison,
						                                         -1, numObjectives, ndCreator,
																 needStrictlyBetter); // -1 signals multi-objective
			}
		}
		else
			return "Unexpected eaType '" + eaType + "'";

		MyGA ga;
		errut::bool_t r;
		if (!(r = ga.run(*creation, *evolver, *calc, stop, popSize, popSize, popSize*2)))
			return "Error running GA: " + r.getErrorString();

		m_best = evolver->getBestIndividuals();
		creation = std::make_shared<eatk::PopulationReuseCreation>(ga.getPopulations());

		return true;
	}
};
