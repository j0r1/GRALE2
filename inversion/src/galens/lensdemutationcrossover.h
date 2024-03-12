#pragma once

#include "graleconfig.h"
#include <eatk/differentialevolutionevolver.h>
#include <vector>

namespace grale
{

class LensDEMutation : public eatk::DifferentialEvolutionMutation
{
public:
	errut::bool_t check(const eatk::Genome &g) override;
	std::shared_ptr<eatk::Genome> mutate(const std::vector<const eatk::Genome*> &genomes, const std::vector<double> &weights) override;
};

class LensDECrossover : public eatk::DifferentialEvolutionCrossover
{
public:
	LensDECrossover(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, bool allowNeg);

	errut::bool_t check(const eatk::Genome &g) override;
	errut::bool_t crossover(double CR, eatk::Genome &mutantDest, const eatk::Genome &origVector) override;
private:
	std::shared_ptr<eatk::RandomNumberGenerator> m_rng;
	const bool m_allowNegative;
};

} // end namespace
