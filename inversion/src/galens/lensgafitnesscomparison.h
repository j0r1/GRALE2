#pragma once

#include "graleconfig.h"
#include <eatk/genomefitness.h>

namespace grale
{

class LensGAFitnessComparison : public eatk::FitnessComparison
{
public:
	LensGAFitnessComparison();
	~LensGAFitnessComparison();
	
	errut::bool_t check(const eatk::Fitness &f) const override;
	bool isFitterThan(const eatk::Fitness &first, const eatk::Fitness &second, size_t objectiveNumber) const override;
};

}