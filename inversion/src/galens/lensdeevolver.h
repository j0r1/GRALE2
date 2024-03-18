#pragma once

#include "graleconfig.h"
#include <eatk/jadeevolver.h>
#include <eatk/differentialevolutionevolver.h>

namespace grale
{

class LensDEEvolver : public eatk::DifferentialEvolutionEvolver
{
public:
	LensDEEvolver(
		const std::shared_ptr<eatk::RandomNumberGenerator> &rng,
		const std::shared_ptr<eatk::DifferentialEvolutionMutation> &mut,
		double F,
		const std::shared_ptr<eatk::DifferentialEvolutionCrossover> &cross,
		double CR,
		const std::shared_ptr<eatk::FitnessComparison> &fitComp, int objectiveNumber = 0, size_t numObjectives = 1,
		const std::shared_ptr<eatk::NonDominatedSetCreator> &ndCreator = nullptr);

	errut::bool_t check(const std::shared_ptr<eatk::Population> &population) override;
	// We need to override this function to copy the calculated scale factors to the genomes
	errut::bool_t createNewPopulation(size_t generation, std::shared_ptr<eatk::Population> &population, size_t targetPopulationSize) override;
};

class LensJADEEvolver : public eatk::JADEEvolver
{
public:
	LensJADEEvolver(
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
