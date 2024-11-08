#pragma once

#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "deparameters.h"
#include "rndparameters.h"
#include "nsga2parameters.h"
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
#include "lensnsga2evolver.h"
#include "randomnumbergenerator.h"
#include "ea.h"
#include "preferredindividualselector.h"
#include "stop.h"
#include "exchange.h"
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

class InversionCommunicatorBase : public InversionCommunicator
{
public:
	InversionCommunicatorBase() { }
	~InversionCommunicatorBase() { }

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
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator0,
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

		std::shared_ptr<grale::RandomNumberGenerator> rng0 = std::make_shared<grale::RandomNumberGenerator>();
		WriteLineStdout("GAMESSAGESTR:RNG SEED: " + std::to_string(rng0->getSeed()));

#if 0
		std::shared_ptr<eatk::RandomNumberGenerator> rng = std::make_shared<RngWrapper>(rng0);
#else
		std::shared_ptr<eatk::RandomNumberGenerator> rng = rng0;
#endif

		std::shared_ptr<grale::LensGAGenomeCalculator> genomeCalculator = std::dynamic_pointer_cast<grale::LensGAGenomeCalculator>(genomeCalculator0);
		if (!genomeCalculator.get())
			return "Calculator does not seem to be of a type derived from LensGAGenomeCalculator";

		auto comparison = genomeCalculator->getFitnessComparison();
		m_selector = std::make_shared<SubsequentBestIndividualSelector>(
								genomeCalculator->getNumberOfObjectives(),
								comparison);

		bool_t r = runGA_next(rng, comparison, popSize, lensFitnessObjectType, calculatorType, calcFactory,
						genomeCalculator0, factoryParamBytes, allEAParams, allConvParams,
						multiPopParams, allEATypes);

		calculatorCleanup();
		return r;
	}

	virtual errut::bool_t runGA_next(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
						 const std::shared_ptr<eatk::FitnessComparison> &comparison,
						 int popSize, const std::string &lensFitnessObjectType,
						 const std::string &calculatorType,
	                     grale::LensGACalculatorFactory &calcFactory, 
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator0,
						 const std::vector<uint8_t> &factoryParamBytes,
						 const std::vector<std::unique_ptr<grale::EAParameters>> &allEAParams,
						 const std::vector<grale::LensGAConvergenceParameters> &allConvParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::vector<std::string> &allEATypes) = 0;

	virtual errut::bool_t getCalculator(const std::string &lensFitnessObjectType,
									const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									eatk::IndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) = 0;

	virtual void calculatorCleanup() { }
	
	std::vector<std::shared_ptr<eatk::Individual>> m_best;
	std::shared_ptr<SubsequentBestIndividualSelector> m_selector;
};
