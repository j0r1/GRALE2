#pragma once

#include "graleconfig.h"
#include "lensinversiongafactorycommon.h"
#include "randomnumbergenerator.h"
#include "vector2d.h"
#include "lensinversionparametersmultiplanegpu.h"
#include "lensinversionbasislensinfo.h"
#include "plummerlensinfo.h"
#include <vector>
#include <memory>

namespace grale
{

class LensInversionParametersMultiPlaneGPU;
class LensInversionGenome;
class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;
class ConfigurationParameters;

// NOTE: the virtual inheritance is again very important!
class GRALE_IMPORTEXPORT LensInversionGAFactoryMultiPlaneGPU : public LensInversionGAFactoryCommon
{
public:
	LensInversionGAFactoryMultiPlaneGPU(std::unique_ptr<LensFitnessObject> fitObj);
	~LensInversionGAFactoryMultiPlaneGPU();

	errut::bool_t init(const LensInversionParametersBase &p) override;

	std::unique_ptr<GravitationalLens> createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const override;

	// These a the EATK functions that we'll override to allow an async calculation
	errut::bool_t startNewCalculation(const eatk::Genome &genome) override;
	errut::bool_t pollCalculate(const eatk::Genome &genome, eatk::Fitness &fitness) override;

	errut::bool_t initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues) override;
	errut::bool_t calculateMassScaleFitness(float scaleFactor, float &fitness) override;
	errut::bool_t calculateTotalFitness(float scaleFactor, float *pFitnessValues) override;
private:
	errut::bool_t analyzeLensBasisFunctions(const std::vector<double> redshifts,
								   const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &basisLenses);
	errut::bool_t analyzeSourceImages(const std::vector<std::shared_ptr<ImagesDataExtended>> &sourceImages,
							 const Cosmology &cosmology,
							 std::vector<std::shared_ptr<ImagesDataExtended>> &imagesVector);
	void convertGenomeSheetValuesToDensities(const std::vector<float> &sheetValues,
											 std::vector<float> &sheetDensities) const;
	void scaleWeights(float scaleFactor);

	std::unique_ptr<LensInversionParametersMultiPlaneGPU> m_currentParams;

	std::vector<float> m_lensRedshifts;
	std::vector<std::vector<PlummerLensInfo>> m_basisLenses;
	std::vector<double> m_basisFunctionMasses;

	std::vector<float> m_sheetDensities;
	std::vector<float> m_sheetMultipliers;
	std::vector<std::vector<float>> m_basePlaneWeights;
	std::vector<std::vector<float>> m_scaledPlaneWeights;

	std::vector<std::shared_ptr<ImagesDataExtended>> m_images;
	std::vector<ImagesDataExtended *> m_reducedImages, m_shortImages;
};

} // end namespace
