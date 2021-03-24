#pragma once

#include "graleconfig.h"
#include "lensgagenomecrossover.h"
#include <mogal2/crossovermutation.h>
#include <mogal2/randomnumbergenerator.h>
#include <cmath>

namespace grale
{

class LensGAIndividual;

class LensGACrossoverBase : public mogal2::PopulationCrossover
{
public:
	LensGACrossoverBase(double beta, bool elitism, bool includeBest, double crossoverRate,
				const std::shared_ptr<mogal2::RandomNumberGenerator> &rng,
				bool allowNegative,
				const std::shared_ptr<mogal2::GenomeMutation> &mutation);

	errut::bool_t check(const std::shared_ptr<mogal2::Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize) override;
protected:
	size_t pickBetaDistIndex(size_t s)
	{
		double x = m_rng->getRandomDouble();
        double val = (1.0-pow(x, 1.0/(1.0+m_beta)))*((double)s);
        size_t r = (size_t)val;
        if (r >= s)
            r = s - 1;
		return r;
	}

	void pickParents(const std::shared_ptr<mogal2::Population> &population, LensGAIndividual **pParent1, LensGAIndividual **pParent2);

	void copyScaleFactorFromFitnessToGenome(const std::shared_ptr<mogal2::Population> &population);
	void copyPopulationIndex(const std::shared_ptr<mogal2::Population> &population);

	errut::bool_t crossover(size_t generation, std::shared_ptr<mogal2::Population> &population, std::shared_ptr<mogal2::Population> &newPop);
	errut::bool_t mutation(size_t mutOffset, std::shared_ptr<mogal2::Population> &newPop);

	virtual LensGAIndividual *pickParent(const std::shared_ptr<mogal2::Population> &population) = 0;
	virtual size_t elitism(std::shared_ptr<mogal2::Population> &population, std::shared_ptr<mogal2::Population> &newPop) = 0;
	
	virtual errut::bool_t sortCheck(const std::shared_ptr<mogal2::Population> &population) = 0;
	virtual errut::bool_t sort(std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize) = 0;

protected:
	std::shared_ptr<mogal2::RandomNumberGenerator> m_rng;
	grale::LensGAGenomeCrossover m_cross;
	double m_beta, m_bestWithMutation, m_bestWithoutMutation, m_crossoverRate;
	std::shared_ptr<mogal2::GenomeMutation> m_mutation;
};

}