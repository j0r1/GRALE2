#pragma once

#include "graleconfig.h"
#include "lensgagenomecalculator.h"
#include "lensfitnessobject.h"
#include <cassert>

namespace grale
{

class LensGAParametricSinglePlaneCalculator : public LensGAGenomeCalculator
{
public:
	LensGAParametricSinglePlaneCalculator(std::unique_ptr<LensFitnessObject> fitObj);
	~LensGAParametricSinglePlaneCalculator();

	errut::bool_t init(const LensInversionParametersBase &params) override;
	errut::bool_t createLens(const eatk::Genome &genome, std::unique_ptr<GravitationalLens> &lens) const override;
	size_t getNumberOfObjectives() const override { assert(m_fitObj.get()); return m_fitObj->getNumberOfFitnessComponents(); }
private:
	std::unique_ptr<LensFitnessObject> m_fitObj;
};

} // end namespace