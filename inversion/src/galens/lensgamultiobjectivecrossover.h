#pragma once

#include "graleconfig.h"
#include "lensgacrossoverbase.h"
#include <eatk/ndsortedpopulation.h>

namespace grale
{

class LensGAMultiObjectiveCrossover : public LensGACrossoverBase
{
public:
	LensGAMultiObjectiveCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
				const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
				bool allowNegative,
				const std::shared_ptr<eatk::GenomeMutation> &mutation,
                size_t numObjectives);

    const std::vector<std::shared_ptr<eatk::Individual>> &getBestIndividuals() const override { return m_sortedPop.getBestIndividuals(); }
private:
	LensGAIndividual *pickParent(const std::shared_ptr<eatk::Population> &population) override;
    void pickParentsRaw(const std::shared_ptr<eatk::Population> &population, LensGAIndividual **pParent1, LensGAIndividual **pParent2) override;
	size_t elitism(std::shared_ptr<eatk::Population> &population, std::shared_ptr<eatk::Population> &newPop) override;
	errut::bool_t sortCheck(const std::shared_ptr<eatk::Population> &population) override;
	errut::bool_t sort(std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize) override;

    eatk::NDSortedPopulation m_sortedPop;
};

}
