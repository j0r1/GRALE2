#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include "gravitationallens.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
#include "retraceparameters.h"
#include <vector>
#include <cassert>

namespace grale
{

class LensInversionParametersParametricSinglePlane : public LensInversionParametersBase
{
public:
	LensInversionParametersParametricSinglePlane();
	LensInversionParametersParametricSinglePlane(
		const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
		double Dd, double zd,
		const GravitationalLens &templateLens,
		double deflectionScale, double potentialScale,
		const std::vector<int> &offsets,
		const std::vector<float> &initMin, const std::vector<float> &initMax,
		const std::vector<float> &hardMin, const std::vector<float> &hardMax,
		bool infOnBoundsViolation,
		const ConfigurationParameters &fitnessObjectParams,
		int devIdx,
		bool randomizeImagePositions,
		uint64_t initialUncertSeed,
		const std::vector<std::pair<size_t, std::string>> &originParameterMapping,
		size_t numOriginParameters,
		bool allowUnusedPriors,
		const std::vector<bool> &retraceImages,
		const std::shared_ptr<TraceParameters> &traceParameters,
		double sourcePlaneDistThreshold, // threshold to accept retrace convergence
		const std::string &clPriorCode,
		bool allowEqualInitRange,
		const std::vector<std::vector<float>> &genomesToCalculateFitnessFor,
		BetaReductionWeightType betaRedWt
		);
	~LensInversionParametersParametricSinglePlane();

	const std::vector<std::shared_ptr<ImagesDataExtended>> &getImages() const { return m_images; }
	double getDd() const { return m_Dd; }
	double getZd() const { return m_zd; }
	double getDeflectionScale() const { return m_deflScale; }
	double getPotentialScale() const { return m_potScale; }
	const GravitationalLens &getTemplateLens() const { assert(m_templateLens.get()); return *m_templateLens; }
	const std::vector<int> &getOffsets() const { return m_offsets; }
	const std::vector<float> &getInitMin() const { return m_initMin; }
	const std::vector<float> &getInitMax() const { return m_initMax; }
	const std::vector<float> &getHardMin() const { return m_hardMin; }
	const std::vector<float> &getHardMax() const { return m_hardMax; }
	bool infinityOnBoundsViolation() const { return m_infOnBoundsViolation; }
	const ConfigurationParameters &getFitnessObjectParameters() const { assert(m_fitObjParams.get()); return *m_fitObjParams; }
	int getDeviceIndex() const { return m_devIdx; }
	bool getRandomizeInputPositions() const { return m_randomizeInputPosition; }
	uint64_t getInitialPositionUncertaintySeed() const { return m_initialUncertSeed; }
	const std::vector<std::pair<size_t, std::string>> &getOriginParameterMapping() const { return m_originParams; }
	size_t getNumberOfOriginParameters() const { return m_numOriginParams; }
	bool shouldAllowUnusedPriors() const { return m_allowUnusedPriors; }
	const std::vector<bool> &shouldRetraceImages() const { return m_retraceImages; }
	const TraceParameters &getRetraceParameters() const { return *m_retraceParams; }
	double getSourcePlaneDistanceThreshold() const { return m_sourcePlaneDistThreshold; }
	std::string getOpenCLPriorCode() const { return m_clPriorCode; }
	bool allowEqualValuesInInitialRange() const { return m_allowEqualInitRange; }
	const std::vector<std::vector<float>> &getGenomesToCalculateFitness() const { return m_genomesToCalculate; }
	BetaReductionWeightType getBetaReductionWeightType() const { return m_betaRedWeigthType; }

	bool write(serut::SerializationInterface &si) const override;
	bool read(serut::SerializationInterface &si) override;
private:
	std::vector<std::shared_ptr<ImagesDataExtended>> m_images;
	double m_Dd = 0, m_zd = 0;
	double m_deflScale = 0, m_potScale = 0;
	std::unique_ptr<GravitationalLens> m_templateLens;
	std::vector<int> m_offsets;
	std::vector<float> m_hardMin, m_hardMax, m_initMin, m_initMax;
	std::unique_ptr<ConfigurationParameters> m_fitObjParams;
	bool m_infOnBoundsViolation = false;
	int m_devIdx = -1;
	bool m_randomizeInputPosition = false;
	uint64_t m_initialUncertSeed = 0;
	std::vector<std::pair<size_t, std::string>> m_originParams;
	size_t m_numOriginParams = 0;
	bool m_allowUnusedPriors = false;
	std::vector<bool> m_retraceImages;
	std::shared_ptr<TraceParameters> m_retraceParams;
	double m_sourcePlaneDistThreshold = 0;
	std::string m_clPriorCode;
	bool m_allowEqualInitRange = false;
	BetaReductionWeightType m_betaRedWeigthType = EqualWeights;

	// Perhaps a better name should be chosen than genome?
	// It's multiple sets of parameters, for which the probability/fitness
	// should be calculated
	std::vector<std::vector<float>> m_genomesToCalculate;
};

} // end namespace
