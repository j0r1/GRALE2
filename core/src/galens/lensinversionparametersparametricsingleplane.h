#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include "gravitationallens.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"
#include <vector>

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
		const ConfigurationParameters &fitnessObjectParams,
		bool uploadFullParameters, int devIdx);
	~LensInversionParametersParametricSinglePlane();

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
	bool m_uploadFillParams = true;
	int m_devIdx = -1;
};

} // end namespace