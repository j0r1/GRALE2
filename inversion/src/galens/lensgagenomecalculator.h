#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include "lensgaindividual.h"
#include "gravitationallens.h"
#include "lensfitnessobject.h"
#include <eatk/genomefitness.h>

namespace grale
{

class LensGAGenomeCalculator : public eatk::GenomeFitnessCalculation
{
public:
    virtual errut::bool_t init(const LensInversionParametersBase &params) = 0;
    virtual errut::bool_t createLens(const LensGAGenome &genome, std::unique_ptr<GravitationalLens> &lens) const = 0;
    virtual size_t getNumberOfObjectives() const = 0;
    virtual bool allowNegativeValues() const = 0;
    virtual size_t getNumberOfBasisFunctions() const = 0;
    virtual size_t getNumberOfSheets() const = 0;
	// TODO: move this to GA parameters
    virtual size_t getMaximumNumberOfGenerations() const = 0;
    // TODO: should not be needed anymore after stop criterion parameters are
    //       moved from lensfitnessgeneral
    virtual const LensFitnessObject &getLensFitnessObject() const = 0;
};

}