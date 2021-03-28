#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include <mogal2/genomefitness.h>

namespace grale
{

class LensGAGenomeCalculator : public mogal2::GenomeFitnessCalculation
{
public:
    virtual errut::bool_t init(const LensInversionParametersBase &params) = 0;
};

}