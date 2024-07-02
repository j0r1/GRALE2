#pragma once

#include "graleconfig.h"
#include <eatk/nsga2evolver.h>

namespace grale
{

class LensNSGA2Evolver : public eatk::NSGA2Evolver
{
public:
	LensNSGA2Evolver(const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
		const std::shared_ptr<eatk::GenomeCrossover> &genomeCrossover,
		const std::shared_ptr<eatk::GenomeMutation> &genomeMutation,
		const std::shared_ptr<eatk::FitnessComparison> &fitComp, size_t numObjectives);
				
	errut::bool_t check(const std::shared_ptr<eatk::Population> &population) override;
	// We need to override this function to copy the calculated scale factors to the genomes
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize) override;
};

}