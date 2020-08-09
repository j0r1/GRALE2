#pragma once

#include "graleconfig.h"
#include "lensinversiongafactorycommon.h"
#include "randomnumbergenerator.h"
#include "mpcudabackprojector.h"
#include "vector2d.h"
#include "lensinversionparametersmultiplanegpu.h"
#include "lensinversionbasislensinfo.h"
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

class LensInversionGAFactoryParamsMultiPlaneGPU : public mogal::GAFactoryParams
{
public:
	LensInversionGAFactoryParamsMultiPlaneGPU();
	LensInversionGAFactoryParamsMultiPlaneGPU(const LensInversionGAFactoryParamsMultiPlaneGPU &src);
	LensInversionGAFactoryParamsMultiPlaneGPU(const LensInversionParametersMultiPlaneGPU &params);
	~LensInversionGAFactoryParamsMultiPlaneGPU();

	const Cosmology &getCosmology() const													{ return m_params.getCosmology(); }
	const std::vector<double> &getLensRedshifts() const										{ return m_params.getLensRedshifts(); }
	const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &getBasisLenses() const { return m_params.getBasisLenses(); }
	const std::vector<std::shared_ptr<ImagesDataExtended>> &getSourceImages() const			{ return m_params.getSourceImages(); }
	double getMassEstimate() const															{ return m_params.getMassEstimate(); }
	bool useMassSheetBasisFunctions() const													{ return m_params.useMassSheetBasisFunctions(); }
	const ConfigurationParameters *getFitnessObjectParameters() const						{ return m_params.getFitnessObjectParameters(); }
	int getMaximumNumberOfGenerations() const												{ return m_params.getMaximumNumberOfGenerations(); }
	bool getAllowNegativeWeights() const													{ return m_params.getAllowNegativeWeights(); }
	const ScaleSearchParameters &getMassScaleSearchParameters() const			   		{ return m_params.getMassScaleSearchParameters(); }
	int getDeviceIndex() const																{ return m_params.getDeviceIndex(); }

	bool write(serut::SerializationInterface &si) const override;
	bool read(serut::SerializationInterface &si) override;
private:
	LensInversionParametersMultiPlaneGPU m_params;
};

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

	bool initializeNewCalculation(const std::vector<float> &basisFunctionWeights, const std::vector<float> &sheetValues) override;
	bool calculateMassScaleFitness(float scaleFactor, float &fitness) override;
	bool calculateTotalFitness(float scaleFactor, float *pFitnessValues) override;
private:
	bool analyzeLensBasisFunctions(const std::vector<double> redshifts,
								   const std::vector<std::vector<std::shared_ptr<LensInversionBasisLensInfo>>> &basisLenses);
	bool analyzeSourceImages(const std::vector<std::shared_ptr<ImagesDataExtended>> &sourceImages,
							 const Cosmology &cosmology,
							 std::vector<std::shared_ptr<ImagesDataExtended>> &imagesVector);
	void convertGenomeSheetValuesToDensities(const std::vector<float> &sheetValues,
											 std::vector<float> &sheetDensities) const;
	bool scaleWeights(float scaleFactor);
	bool checkCUDAInit();

	std::unique_ptr<LensInversionGAFactoryParamsMultiPlaneGPU> m_currentParams;
	std::shared_ptr<MPCUDABackProjector> m_cudaBpShort, m_cudaBpFull;

	std::vector<float> m_lensRedshifts;
	std::vector<std::vector<PlummerLensInfo>> m_basisLenses;
	std::vector<double> m_basisFunctionMasses;

	std::vector<float> m_sheetDensities;
	std::vector<float> m_sheetMultipliers;
	std::vector<std::vector<float>> m_basePlaneWeights;
	std::vector<std::vector<float>> m_scaledPlaneWeights;

	std::vector<std::shared_ptr<ImagesDataExtended>> m_images;
	std::vector<ImagesDataExtended *> m_reducedImages, m_shortImages;
	std::string m_libraryPath;
	bool m_cudaInitialized;
	bool m_cudaInitAttempted;
};

} // end namespace
