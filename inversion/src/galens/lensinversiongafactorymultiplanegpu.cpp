#include "lensinversiongafactorymultiplanegpu.h"

namespace grale
{

LensInversionGAFactoryMultiPlaneGPU::LensInversionGAFactoryMultiPlaneGPU()
{

}

LensInversionGAFactoryMultiPlaneGPU::~LensInversionGAFactoryMultiPlaneGPU()
{

}

mogal::GAFactoryParams *LensInversionGAFactoryMultiPlaneGPU::createParamsInstance() const
{
    return nullptr;
}

const mogal::GAFactoryParams *LensInversionGAFactoryMultiPlaneGPU::getCurrentParameters() const
{
    return nullptr;
}

bool LensInversionGAFactoryMultiPlaneGPU::init(const mogal::GAFactoryParams *p)
{
    return true;
}

GravitationalLens *LensInversionGAFactoryMultiPlaneGPU::createLens(const LensInversionGenome &genome, std::string &errStr) const
{
    return nullptr;
}

bool LensInversionGAFactoryMultiPlaneGPU::initializeNewCalculation(const std::vector<float> &masses, const std::vector<float> &sheetValues)
{
    return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::calculateMassScaleFitness(float scaleFactor, float &fitness)
{
    return true;
}

bool LensInversionGAFactoryMultiPlaneGPU::calculateTotalFitness(float scaleFactor, float *pFitnessValues)
{
    return true;
}

} // enc namespace