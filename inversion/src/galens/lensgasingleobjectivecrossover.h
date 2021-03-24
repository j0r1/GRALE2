#pragma once

#include "graleconfig.h"
#include "lensgacrossoverbase.h"
#include <mogal2/simplesortedpopulation.h>

namespace grale
{

class LensGASingleObjectiveCrossover : public LensGACrossoverBase
{
public:
	LensGASingleObjectiveCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
				const std::shared_ptr<mogal2::RandomNumberGenerator> &rng,
				bool allowNegative,
				const std::shared_ptr<mogal2::GenomeMutation> &mutation);
    
    const std::vector<std::shared_ptr<mogal2::Individual>> &getBestIndividuals() const override { return m_sortedPop.getBestIndividuals(); }
private:
	LensGAIndividual *pickParent(const std::shared_ptr<mogal2::Population> &population) override;
	size_t elitism(std::shared_ptr<mogal2::Population> &population, std::shared_ptr<mogal2::Population> &newPop) override;
	errut::bool_t sortCheck(const std::shared_ptr<mogal2::Population> &population) override;
	errut::bool_t sort(std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize) override;

	mogal2::SimpleSortedPopulation m_sortedPop;
};

}
