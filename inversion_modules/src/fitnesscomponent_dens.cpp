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
	addRecognizedTypeName("kappathresholdpoints");
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
	if (typeName != "kappathresholdpoints")
		return true; // not relevant, ignore

	string thresStr;
	double threshold = 1.0;

	if (imgDat.hasExtraParameter("threshold"))
	{
		if (!imgDat.getExtraParameter("threshold", threshold))
		{
			setErrorString("Extra parameter 'threshold' is present, but doesn't seem to be a real value");
			return false;
		}
	}

	if (threshold < 0)
	{
		setErrorString("The specified threshold is negative");
		return false;
	}

	if (imgDat.getNumberOfImages() != 1)
	{
		setErrorString("Only a single large 'image' (containing a set of points for which kappa should be inspected) is allowed");
		return false;
	}
	if (imgDat.getNumberOfImagePoints(0) < 1)
	{
		setErrorString("No image points are present");
		return false;
	}

	needCalcDeflDeriv = true;
	needCalcConvergence = true;

	addImagesDataIndex(idx);
	m_thresholds.push_back((float)threshold);
	return true;
}

bool FitnessComponent_KappaThreshold::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateKappaThresholdFitness(iface, getUsedImagesDataIndices(), m_thresholds);
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
