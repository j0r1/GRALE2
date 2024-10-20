#include "lensinversionparametersparametricsingleplane.h"
#include <iostream>

using namespace std;

namespace grale
{

LensInversionParametersParametricSinglePlane::LensInversionParametersParametricSinglePlane()
{
}

LensInversionParametersParametricSinglePlane::LensInversionParametersParametricSinglePlane(
		const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
		double Dd, double zd,
		const GravitationalLens &templateLens,
		double deflectionScale, double potentialScale,
		const std::vector<int> &offsets,
		const std::vector<float> &initMin, const std::vector<float> &initMax,
		const std::vector<float> &hardMin, const std::vector<float> &hardMax,
		const ConfigurationParameters &fitnessObjectParams,
		bool uploadFullParameters, int devIdx)
{
	m_images = images;
	m_Dd = Dd;
	m_zd = zd;
	m_templateLens = move(templateLens.createCopy());
	m_deflScale = deflectionScale;
	m_potScale = potentialScale;
	m_offsets = offsets;
	m_initMin = initMin;
	m_initMax = initMax;
	m_hardMin = hardMin;
	m_hardMax = hardMax;
	m_fitObjParams = make_unique<ConfigurationParameters>(fitnessObjectParams);
	m_uploadFillParams = uploadFullParameters;
	m_devIdx = devIdx;
}

LensInversionParametersParametricSinglePlane::~LensInversionParametersParametricSinglePlane()
{
}

bool LensInversionParametersParametricSinglePlane::write(serut::SerializationInterface &si) const
{
	if (m_offsets.size() != m_initMin.size() || m_offsets.size() != m_initMax.size() ||
		m_offsets.size() != m_hardMin.size() || m_offsets.size() != m_hardMax.size())
	{
		setErrorString("Incompatible lengths of offsets, initMin, initMax, hardMin, hardMax");
		return false;
	}

	if (!m_templateLens || !m_fitObjParams)
	{
		setErrorString("Lens or fitness object parameters are not set");
		return false;
	}

	int32_t numImgs = m_images.size();
	if (!si.writeInt32(numImgs))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	for (auto &img : m_images)
	{
		if (!img->write(si))
		{
			setErrorString("Can't write image: " + img->getErrorString());
			return false;
		}
	}

	double dParams[] = { m_Dd, m_zd, m_deflScale, m_potScale };
	if (!si.writeDoubles(dParams, 4))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!m_templateLens->write(si))
	{
		setErrorString("Can't write template lens: " + m_templateLens->getErrorString());
		return false;
	}

	int32_t numOffsets = m_offsets.size();

	vector<int32_t> offsets;
	for (auto x : m_offsets)
		offsets.push_back((int32_t)x);

	if (!si.writeInt32(numOffsets) ||
		!si.writeInt32s(offsets) ||
		!si.writeFloats(m_initMin) ||
		!si.writeFloats(m_initMax) ||
		!si.writeFloats(m_hardMin) ||
		!si.writeFloats(m_hardMax)
		)
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (!m_fitObjParams->write(si))
	{
		setErrorString(m_fitObjParams->getErrorString());
		return false;
	}

	int32_t iParams[] = { (m_uploadFillParams)?1:0, (int32_t)m_devIdx };
	if (!si.writeInt32s(iParams, 2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	return true;
}

bool LensInversionParametersParametricSinglePlane::read(serut::SerializationInterface &si)
{
	int32_t numImgs;
	if (!si.readInt32(&numImgs))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (numImgs < 0)
	{
		setErrorString("Read invalid number of images");
		return false;
	}

	m_images.clear();
	for (int32_t i = 0 ; i < numImgs ; i++)
	{
		shared_ptr<ImagesDataExtended> img = make_shared<ImagesDataExtended>();
		if (!img->read(si))
		{
			setErrorString("Can't read image: " + img->getErrorString());
			return false;
		}
		m_images.push_back(img);
	}
	double dParams[4];
	if (!si.readDoubles(dParams, 4))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	m_Dd = dParams[0];
	m_zd = dParams[1];
	m_deflScale = dParams[2];
	m_potScale = dParams[3];

	string errStr;
	if (!GravitationalLens::read(si, m_templateLens, errStr))
	{
		setErrorString("Can't read template lens: " + errStr);
		return false;
	}

	int32_t numOffsets;
	if (!si.readInt32(&numOffsets))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (numOffsets < 0)
	{
		setErrorString("Read invalid size of offsets vector");
		return false;
	}

	vector<int32_t> offsets(numOffsets);
	m_initMin.resize(numOffsets);
	m_initMax.resize(numOffsets);
	m_hardMin.resize(numOffsets);
	m_hardMax.resize(numOffsets);

	if (!si.readInt32s(offsets) ||
		!si.readFloats(m_initMin) ||
		!si.readFloats(m_initMax) ||
		!si.readFloats(m_hardMin) ||
		!si.readFloats(m_hardMax)
	)
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_offsets.clear();
	for (int32_t x : offsets)
		m_offsets.push_back((int)x);

	m_fitObjParams = make_unique<ConfigurationParameters>();
	if (!m_fitObjParams->read(si))
	{
		setErrorString(m_fitObjParams->getErrorString());
		return false;
	}

	int32_t iParams[2];
	if (!si.readInt32s(iParams, 2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_uploadFillParams = (iParams[0] == 0)?false:true;
	m_devIdx = (int)iParams[1];

	return true;
}

}
