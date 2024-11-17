#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include "gravitationallens.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
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
		bool uploadFullParameters, int devIdx,
		bool randomizeImagePositions,
		uint64_t initialUncertSeed
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
	bool alwaysUploadFullParameters() const { return m_uploadFillParams; }
	int getDeviceIndex() const { return m_devIdx; }
	bool getRandomizeInputPositions() const { return m_randomizeInputPosition; }
	uint64_t getInitialPositionUncertaintySeed() const { return m_initialUncertSeed; }

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
	bool m_uploadFillParams = true;
	int m_devIdx = -1;
	bool m_randomizeInputPosition = false;
	uint64_t m_initialUncertSeed = 0;
};

} // end namespace
