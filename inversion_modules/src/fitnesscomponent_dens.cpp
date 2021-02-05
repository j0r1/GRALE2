#include "fitnesscomponent_dens.h"
#include "fitnessutil.h"
#include <grale/constants.h>
#include <grale/imagesdataextended.h>
#include <grale/projectedimagesinterface.h>

using namespace std;
using namespace errut;

namespace grale
{

// FitnessComponent_KappaThreshold

FitnessComponent_KappaThreshold::FitnessComponent_KappaThreshold(FitnessComponentCache *pCache) 
	: FitnessComponent("kappathreshold", pCache)
{
	addRecognizedTypeName("kappathresholdmin");
	addRecognizedTypeName("kappathresholdmax");
}

FitnessComponent_KappaThreshold::~FitnessComponent_KappaThreshold()
{
}

bool FitnessComponent_KappaThreshold::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "kappathresholdmin" && typeName != "kappathresholdmax")
		return true; // not relevant, ignore

	if (!imgDat.hasProperty(ImagesData::Kappa))
	{
		setErrorString("image data should have 'kappa' property");
		return false;
	}

	needCalcDeflDeriv = true;
	needCalcConvergence = true;

	addImagesDataIndex(idx);
	m_minMax.push_back((typeName == "kappathresholdmin")?false:true);
	return true;
}

bool FitnessComponent_KappaThreshold::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	const vector<int> &srcs = getUsedImagesDataIndices();
	assert(srcs.size() == m_minMax.size());

	float diffSum = 0;

	for (size_t sIdx = 0 ; sIdx < srcs.size() ; sIdx++)
	{
		bool isMax = m_minMax[sIdx];
		const float isMaxFactor = (isMax)?1.0f:-1.0f;

		const int s = srcs[sIdx];
		assert(s >= 0 && s < iface.getNumberOfSources());

		const int numPoints = iface.getNumberOfImagePoints(s);
		const float *pKappaOrig = iface.getOriginalProperties(ImagesData::Kappa, s);
		const float *pKappa = iface.getConvergence(s);

		for (int i = 0 ; i < numPoints ; i++)
		{
			float diff = (pKappa[i] - pKappaOrig[i])*isMaxFactor;
			if (diff > 0)
				diffSum += diff;
		}
	}

	fitness = diffSum;
	return true;
}

// FitnessComponent_KappaGradient

FitnessComponent_KappaGradient::FitnessComponent_KappaGradient(FitnessComponentCache *pCache)
	: FitnessComponent("kappagradient", pCache),
	  m_srcIdx(-1)
{
	addRecognizedTypeName("kappagradientgrids");
}

FitnessComponent_KappaGradient::~FitnessComponent_KappaGradient()
{
}

bool FitnessComponent_KappaGradient::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
								bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
								bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "kappagradientgrids") // ignore
		return true;

	if (m_srcIdx >= 0)
	{
		setErrorString("A grid for the kappa gradient has already been set, only one can be present");
		return false;
	}

	bool_t r;
	if (!(r = m_gradCalc.check(imgDat)))
	{
		setErrorString(r.getErrorString());
		return false;
	}
	
	needCalcDeflDeriv = true; // needed for convergence calculation
	needCalcConvergence = true;

	m_srcIdx = idx;
	addImagesDataIndex(idx);

	return true;
}

bool FitnessComponent_KappaGradient::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	const vector<float> &gradientSizes = m_gradCalc.getGradientSquaredSizes(m_srcIdx, iface);
	float sum = 0;

	for (auto x : gradientSizes)
		sum += x;
	
	float factor = (float)(iface.getAngularScale()/ANGLE_ARCSEC);

	fitness = (sum/gradientSizes.size())/(factor*factor);
	return true;
}

} // end namespace
