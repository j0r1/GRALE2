#pragma once

#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactorycommon.h"
#include "inversioncommunicatornewga.h"
#include "gslrngwrapper.h"
#include "constants.h"
#include "lensgaindividual.h"
#include "lensgagenomemutation.h"
#include "lensgagenomecrossover.h"
#include "lensgafitnesscomparison.h"
#include "lensgasingleobjectivecrossover.h"
#include "lensgamultiobjectivecrossover.h"
#include "lensgastopcriterion.h"
#include <serut/memoryserializer.h>
#include <eatk/evolutionaryalgorithm.h>
#include <eatk/singlethreadedpopulationfitnesscalculation.h>
#include <eatk/multithreadedpopulationfitnesscalculation.h>
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

	errut::bool_t onBeforeFitnessCalculation(size_t generation, std::shared_ptr<eatk::Population> &population)
	{
		m_timer.start();
	// 	cout << "# Generation " << generation << ", before calculation: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, std::shared_ptr<eatk::Population> &population)
	{
		m_timer.stop();
		m_intervals.push_back(m_timer.duration());
	// 	cout << "# Generation " << generation << ", calculated: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
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
	Stop(size_t maxGenerations, const std::shared_ptr<grale::LensGAGenomeMutation> &mutation)
		: grale::LensGAStopCriterion(maxGenerations, mutation) { }
protected:
	void onReport(const std::string &s)	const override
	{
		WriteLineStdout("GAMESSAGESTR:" + s);
	}
};

// TODO: rename this, is from copy-paste
class NewGACommunicatorBase : public InversionCommunicator
{
public:
	NewGACommunicatorBase() { }
	~NewGACommunicatorBase() { }
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
						 const grale::GAParameters &params)
	{

		errut::bool_t r;

		m_selector = std::make_shared<SubsequentBestIndividualSelector>(
								genomeCalculator->getNumberOfObjectives(),
								std::make_shared<grale::LensGAFitnessComparison>());

		std::shared_ptr<GslRNGWrapper> rng = std::make_shared<GslRNGWrapper>();
		MyGA ga;

		WriteLineStdout("GAMESSAGESTR:RNG SEED: " + std::to_string(rng->getSeed()));

		auto mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 
					   1.0, // chance multiplier; has always been set to one
					   genomeCalculator->allowNegativeValues(),
					   0, // mutation amplitude, will be in the stop criterion
					   true); // absolute or small mutation, will be set in the stop criterion

		grale::LensGAIndividualCreation creation(rng, 
						  genomeCalculator->getNumberOfBasisFunctions(),
						  genomeCalculator->getNumberOfSheets(),
		                  genomeCalculator->allowNegativeValues(),
						  genomeCalculator->getNumberOfObjectives());

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

		std::shared_ptr<eatk::PopulationFitnessCalculation> calc;
		if (!(r = getCalculator(lensFitnessObjectType, calculatorType, calcFactory, genomeCalculator,
								factoryParamBytes, creation, calc)))
			return "Can't get calculator: " + r.getErrorString();
		
		Stop stop(genomeCalculator->getMaximumNumberOfGenerations(), mutation);
		if (!(r = stop.initialize(genomeCalculator->getLensFitnessObject())))
		{
			calculatorCleanup();
			return "Error initializing convergence checker: " + r.getErrorString();
		}

		if (!(r = ga.run(creation, *cross, *calc, stop, popSize)))
		{
			calculatorCleanup();
			return "Error running GA: " + r.getErrorString();
		}

		calculatorCleanup();

		m_best = cross->getBestIndividuals();

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
