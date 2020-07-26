#pragma once

#include "graleconfig.h"
#include "lensinversiongafactorycommon.h"
#include "randomnumbergenerator.h"
#include "backprojectmatrixnew.h"
#include "vector2d.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include <vector>
#include <memory>

namespace grale
{

class LensInversionGAFactoryParamsMultiPlaneGPU;
class LensInversionGenome;
class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;
class ConfigurationParameters;

// NOTE: the virtual inheritance is again very important!
class GRALE_IMPORTEXPORT LensInversionGAFactoryMultiPlaneGPU : public virtual LensInversionGAFactoryCommon
{
public:
	LensInversionGAFactoryMultiPlaneGPU();
	~LensInversionGAFactoryMultiPlaneGPU();

	mogal::GAFactoryParams *createParamsInstance() const override;
	const mogal::GAFactoryParams *getCurrentParameters() const override;

	bool init(const mogal::GAFactoryParams *p) override;

	GravitationalLens *createLens(const LensInversionGenome &genome, std::string &errStr) const override;

	bool initializeNewCalculation(const std::vector<float> &masses, const std::vector<float> &sheetValues) override;
	bool calculateMassScaleFitness(float scaleFactor, float &fitness) override;
	bool calculateTotalFitness(float scaleFactor, float *pFitnessValues) override;
private:
};

} // end namespace
