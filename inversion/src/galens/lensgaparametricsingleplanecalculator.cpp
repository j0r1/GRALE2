#include "lensgaparametricsingleplanecalculator.h"

using namespace errut;
using namespace std;

namespace grale
{

LensGAParametricSinglePlaneCalculator::LensGAParametricSinglePlaneCalculator(std::unique_ptr<LensFitnessObject> fitObj)
	: m_fitObj(move(fitObj))
{
}

LensGAParametricSinglePlaneCalculator::~LensGAParametricSinglePlaneCalculator()
{
}

bool_t LensGAParametricSinglePlaneCalculator::init(const LensInversionParametersBase &params)
{
	return "TODO";
}

bool_t LensGAParametricSinglePlaneCalculator::createLens(const eatk::Genome &genome, std::unique_ptr<GravitationalLens> &lens) const
{
	return "TODO";
}

} // end namespace