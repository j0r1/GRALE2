#pragma once

#include "graleconfig.h"
#include "lensgagenomecrossover.h"
#include <mogal2/crossovermutation.h>
#include <mogal2/randomnumbergenerator.h>
#include <mogal2/simplesortedpopulation.h>

namespace grale
{

class LensGASingleFitnessCrossover : public mogal2::PopulationCrossover
{
public:
	LensGASingleFitnessCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
				const std::shared_ptr<mogal2::RandomNumberGenerator> &rng,
				bool allowNegative,
				const std::shared_ptr<mogal2::GenomeMutation> &mutation);

	errut::bool_t check(const std::shared_ptr<mogal2::Population> &population) override;
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize) override;
    
    const std::vector<std::shared_ptr<mogal2::Individual>> &getBestIndividuals() const override { return m_sortedPop.getBestIndividuals(); }
private:
	void copyScaleFactorFromFitnessToGenome(std::shared_ptr<mogal2::Population> &population);
	size_t elitism(std::shared_ptr<mogal2::Population> &population, std::shared_ptr<mogal2::Population> &newPop);
	errut::bool_t crossover(size_t generation, std::shared_ptr<mogal2::Population> &population, std::shared_ptr<mogal2::Population> &newPop);
	errut::bool_t mutation(size_t mutOffset, std::shared_ptr<mogal2::Population> &newPop);

	std::shared_ptr<mogal2::RandomNumberGenerator> m_rng;
	mogal2::SimpleSortedPopulation m_sortedPop;
	grale::LensGAGenomeCrossover m_cross;
	double m_beta, m_bestWithMutation, m_bestWithoutMutation, m_crossoverRate;
	std::shared_ptr<mogal2::GenomeMutation> m_mutation;
};

}
