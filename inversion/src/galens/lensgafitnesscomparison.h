#pragma once

#include "graleconfig.h"
#include <mogal2/genomefitness.h>

namespace grale
{

class LensGAFitnessComparison : public mogal2::FitnessComparison
{
public:
	LensGAFitnessComparison();
	~LensGAFitnessComparison();
	
	errut::bool_t check(const mogal2::Fitness &f) const override;
	bool isFitterThan(const mogal2::Fitness &first, const mogal2::Fitness &second, size_t objectiveNumber) const override;
};

}