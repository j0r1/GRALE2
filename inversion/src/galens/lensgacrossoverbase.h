#pragma once

#include "graleconfig.h"
#include "lensgagenomecrossover.h"
#include <eatk/crossovermutation.h>
#include <eatk/randomnumbergenerator.h>
#include <cmath>
#include <cassert>

namespace grale
{

class LensGAIndividual;

class LensGACrossoverBase : public eatk::PopulationEvolver
{
public:
	LensGACrossoverBase(double beta, bool elitism, bool includeBest, double crossoverRate,
				const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
				bool allowNegative,
				const std::shared_ptr<eatk::GenomeMutation> &mutation);

	errut::bool_t check(const std::shared_ptr<eatk::Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize) override;
protected:
	size_t pickRandomNumber(size_t s)
	{
		assert(s > 0);
		size_t i = (size_t)(m_rng->getRandomDouble()*(double)s);
		if (i >= s)
			i = s-1;
		assert(i >= 0 && i < s);
		return i;
	}
	size_t pickBetaDistIndex(size_t s)
	{
		double x = m_rng->getRandomDouble();
		double val = (1.0-pow(x, 1.0/(1.0+m_beta)))*((double)s);
		size_t r = (size_t)val;
		if (r >= s)
			r = s - 1;
		return r;
	}

	void pickParentsNoInbreed(const std::shared_ptr<eatk::Population> &population, LensGAIndividual **pParent1, LensGAIndividual **pParent2);

	void copyScaleFactorFromFitnessToGenome(const std::shared_ptr<eatk::Population> &population);
	void copyPopulationIndex(const std::shared_ptr<eatk::Population> &population);

	errut::bool_t crossover(size_t generation, std::shared_ptr<eatk::Population> &population, std::shared_ptr<eatk::Population> &newPop);
	errut::bool_t mutation(size_t mutOffset, std::shared_ptr<eatk::Population> &newPop);

	virtual size_t elitism(std::shared_ptr<eatk::Population> &population, std::shared_ptr<eatk::Population> &newPop) = 0;
	virtual LensGAIndividual *pickParent(const std::shared_ptr<eatk::Population> &population) = 0;
	virtual void pickParentsRaw(const std::shared_ptr<eatk::Population> &population, LensGAIndividual **pParent1, LensGAIndividual **pParent2) = 0;
	virtual errut::bool_t sortCheck(const std::shared_ptr<eatk::Population> &population) = 0;
	virtual errut::bool_t sort(std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize) = 0;
protected:
	std::shared_ptr<eatk::RandomNumberGenerator> m_rng;
	grale::LensGAGenomeCrossover m_cross;
	double m_beta, m_bestWithMutation, m_bestWithoutMutation, m_crossoverRate;
	std::shared_ptr<eatk::GenomeMutation> m_mutation;
};

}