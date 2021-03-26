#pragma once

#include "log.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
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
#include "galensmodule.h"
#include "lensgastopcriterion.h"
#include "lensfitnesscalculation.h"
#include <serut/memoryserializer.h>
#include <mogal/geneticalgorithm.h>
#include <mogal2/geneticalgorithm.h>
#include <mogal2/singlethreadedpopulationfitnesscalculation.h>
#include <mogal2/multithreadedpopulationfitnesscalculation.h>
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


class MyGA : public mogal2::GeneticAlgorithm
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

	errut::bool_t onBeforeFitnessCalculation(size_t generation, std::shared_ptr<mogal2::Population> &population)
	{
		m_timer.start();
	// 	cout << "# Generation " << generation << ", before calculation: " << endl;
	// 	for (auto &i : population->individuals())
	// 		cout << i->fitness()->toString() << endl;
		return true;
	}

    errut::bool_t onFitnessCalculated(size_t generation, std::shared_ptr<mogal2::Population> &population)
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
	virtual errut::bool_t select(const std::vector<std::shared_ptr<mogal2::Individual>> &best,
								 std::shared_ptr<mogal2::Individual> &selected)
	{
		return "Not implemented";
	}
};

class SubsequentBestIndividualSelector : public PreferredIndividualSelector
{
public:
	SubsequentBestIndividualSelector(size_t numObjectives,
									 const std::shared_ptr<mogal2::FitnessComparison> &fitComp)
		: m_numObjectives(numObjectives), m_cmp(fitComp) { }
	~SubsequentBestIndividualSelector() { }
	
	errut::bool_t select(const std::vector<std::shared_ptr<mogal2::Individual>> &bestIndividuals,
								 std::shared_ptr<mogal2::Individual> &selected) override
	{
		std::vector<std::shared_ptr<mogal2::Individual>> genomes = bestIndividuals;
		std::vector<std::shared_ptr<mogal2::Individual>> genomes2;
		std::shared_ptr<mogal2::Individual> best;

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
	std::shared_ptr<mogal2::FitnessComparison> m_cmp;
};

// TODO: rename this, is from copy-paste
class NewGACommunicatorBase : public InversionCommunicator
{
public:
	NewGACommunicatorBase() { }
	~NewGACommunicatorBase() { }
protected:	
	bool hasGeneticAlgorithm() const override { return true; }
	void getAllBestGenomes(std::vector<std::shared_ptr<mogal2::Individual>> &bestGenomes) override
	{
		bestGenomes = m_best;
	}

	std::shared_ptr<mogal2::Individual> getPreferredBestGenome() override
	{
		if (m_best.size() == 0)
			return nullptr;
		if (!m_selector.get())
			return nullptr;
		
		errut::bool_t r;
		std::shared_ptr<mogal2::Individual> selected;
		if (!(r = m_selector->select(m_best, selected)))
		{
			std::cerr << "Error in preferred individual selection: " << r.getErrorString() << std::endl;
			return nullptr;
		}

		return selected;
	}

	errut::bool_t runGA(int popSize, mogal::GAFactory &factory, grale::GAParameters &params,
	             const std::string &moduleDir, const std::string &moduleFile, grale::GALensModule &module,
	             const std::vector<uint8_t> &factoryParamBytes) override
	{
		errut::bool_t r;
		grale::LensInversionGAFactoryCommon &gaFactory = dynamic_cast<grale::LensInversionGAFactoryCommon&>(factory);
		size_t numObjectives = gaFactory.getNumberOfFitnessComponents();

		m_selector = std::make_shared<SubsequentBestIndividualSelector>(numObjectives,
								std::make_shared<grale::LensGAFitnessComparison>());

		std::shared_ptr<mogal2::RandomNumberGenerator> rng = std::make_shared<GslRNGWrapper>();
		MyGA ga;

		auto mutation = std::make_shared<grale::LensGAGenomeMutation>(rng, 
					   gaFactory.getChanceMultiplier(),
					   gaFactory.allowNegativeValues(),
					   gaFactory.getMutationAmplitude(),
					   gaFactory.useAbsoluteMutation());

		grale::LensGAIndividualCreation creation(rng, gaFactory.getNumberOfBasisFunctions(), gaFactory.getNumberOfSheets(),
		                  gaFactory.allowNegativeValues(), gaFactory.getNumberOfFitnessComponents());

		std::shared_ptr<grale::LensGACrossoverBase> cross;
		if (numObjectives == 1)
			cross = std::make_shared<grale::LensGASingleObjectiveCrossover>(params.getSelectionPressure(),
					      params.getUseElitism(),
						  params.getAlwaysIncludeBest(),
						  params.getCrossOverRate(),
						  rng,
						  gaFactory.allowNegativeValues(),
						  mutation);
		else
			cross = std::make_shared<grale::LensGAMultiObjectiveCrossover>(params.getSelectionPressure(),
					      params.getUseElitism(),
						  params.getAlwaysIncludeBest(),
						  params.getCrossOverRate(),
						  rng,
						  gaFactory.allowNegativeValues(),
						  mutation,
						  numObjectives);


		std::shared_ptr<mogal2::PopulationFitnessCalculation> calc;
		if (!(r = getCalculator(gaFactory, moduleDir, moduleFile, module, factoryParamBytes, creation, calc)))
			return "Can't get calculator: " + r.getErrorString();
		
		grale::LensGAStropCriterion stop(gaFactory.getMaximumNumberOfGenerations(), mutation);
		if (!(r = stop.initialize(gaFactory.getFitnessObject())))
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
		std::cout << "Best: " << std::endl;
		for (auto &b: m_best)
			std::cout << b->fitness()->toString() << std::endl;

		return true;
	}

	virtual errut::bool_t getCalculator(grale::LensInversionGAFactoryCommon &gaFactory, 
								const std::string &moduleDir, const std::string &moduleFile, grale::GALensModule &module,
	             				const std::vector<uint8_t> &factoryParamBytes,
								grale::LensGAIndividualCreation &creation,
								std::shared_ptr<mogal2::PopulationFitnessCalculation> &calc) = 0;

	virtual void calculatorCleanup() { }
	
	std::vector<std::shared_ptr<mogal2::Individual>> m_best;
	std::shared_ptr<SubsequentBestIndividualSelector> m_selector;
};
