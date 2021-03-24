#pragma once

#include "graleconfig.h"
#include "lensgacrossoverbase.h"
#include <mogal2/ndsortedpopulation.h>

namespace grale
{

class LensGAMultiObjectiveCrossover : public LensGACrossoverBase
{
public:
	LensGAMultiObjectiveCrossover(double beta, bool elitism, bool includeBest, double crossoverRate,
				const std::shared_ptr<mogal2::RandomNumberGenerator> &rng,
				bool allowNegative,
				const std::shared_ptr<mogal2::GenomeMutation> &mutation,
                size_t numObjectives);

    const std::vector<std::shared_ptr<mogal2::Individual>> &getBestIndividuals() const override { return m_sortedPop.getBestIndividuals(); }
private:
	LensGAIndividual *pickParent(const std::shared_ptr<mogal2::Population> &population) override;
    void pickParentsRaw(const std::shared_ptr<mogal2::Population> &population, LensGAIndividual **pParent1, LensGAIndividual **pParent2) override;
	size_t elitism(std::shared_ptr<mogal2::Population> &population, std::shared_ptr<mogal2::Population> &newPop) override;
	errut::bool_t sortCheck(const std::shared_ptr<mogal2::Population> &population) override;
	errut::bool_t sort(std::shared_ptr<mogal2::Population> &population, size_t targetPopulationSize) override;

    mogal2::NDSortedPopulation m_sortedPop;
};

}
