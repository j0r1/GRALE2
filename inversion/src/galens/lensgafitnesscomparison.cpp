#include "lensgafitnesscomparison.h"
#include "lensgaindividual.h"

using namespace errut;

namespace grale
{

LensGAFitnessComparison::LensGAFitnessComparison()
{
}

LensGAFitnessComparison::~LensGAFitnessComparison()
{
}

bool_t LensGAFitnessComparison::check(const eatk::Fitness &f) const
{
	if (!dynamic_cast<const LensGAFitness*>(&f))
		return "Fitness is of wrong type";
	return true;
}

bool LensGAFitnessComparison::isFitterThan(const eatk::Fitness &first, const eatk::Fitness &second, size_t objectiveNumber) const
{
	const LensGAFitness &f1 = static_cast<const LensGAFitness &>(first);
	const LensGAFitness &f2 = static_cast<const LensGAFitness &>(second);

	assert(f1.m_fitnesses.size() == f2.m_fitnesses.size());
	assert(objectiveNumber < f1.m_fitnesses.size());

	return f1.m_fitnesses[objectiveNumber] < f2.m_fitnesses[objectiveNumber];
}

}