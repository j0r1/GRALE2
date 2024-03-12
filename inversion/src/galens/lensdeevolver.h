#pragma once

#include "graleconfig.h"
#include <eatk/jadeevolver.h>

namespace grale
{

class LensDEEvolver : public eatk::JADEEvolver
{
public:
	LensDEEvolver(
		const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
		const std::shared_ptr<eatk::DifferentialEvolutionMutation> &mut,
		const std::shared_ptr<eatk::DifferentialEvolutionCrossover> &cross,
		const std::shared_ptr<eatk::FitnessComparison> &fitComp
		// TODO: other parameters
		);

	errut::bool_t check(const std::shared_ptr<eatk::Population> &population) override;
	// We need to override this function to copy the calculated scale factors to the genomes
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize) override;
};

} // end namespace
