#include "fitnesscomponent_defl.h"
#include "fitnessutil.h"
#include "constants.h"
#include "imagesdataextended.h"
#include "projectedimagesinterface.h"

using namespace std;
using namespace errut;

namespace grale
{

// FitnessComponent_DeflectionAngle

FitnessComponent_DeflectionAngle::FitnessComponent_DeflectionAngle(const std::shared_ptr<FitnessComponentCache> &pCache)
	: FitnessComponent("deflectionangle", pCache)
{
	addRecognizedTypeName("deflectionangles");
}

FitnessComponent_DeflectionAngle::~FitnessComponent_DeflectionAngle()
{
}

bool FitnessComponent_DeflectionAngle::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
								bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
								bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "deflectionangles") // ignore
		return true;

	if (!imgDat.hasProperty(ImagesData::DeflectionComponent1) || !imgDat.hasProperty(ImagesData::DeflectionComponent2))
	{
		setErrorString("The deflection angle property needs to be set");
		return false;
	}

	if (imgDat.getDs() != 1.0 || imgDat.getDds() != 1.0)
	{
		setErrorString("For deflection angle data, Dds and Ds need to be set to 1");
		return false;
	}

	m_srcIdx.push_back(idx);
	addImagesDataIndex(idx);
	needCalcDeflections = true;

	return true;
}

bool FitnessComponent_DeflectionAngle::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	double unitScale = (iface.getAngularScale()/ANGLE_ARCSEC);
	float deflAngleDiffSum2 = 0;
	int totalNumPts = 0;

	for (int s : m_srcIdx)
	{
		assert(s >= 0 && s < iface.getNumberOfSources());
		int numPoints = iface.getNumberOfImagePoints(s);
		totalNumPts += numPoints;

		assert(iface.hasOriginalProperty(ImagesData::DeflectionComponent1, s) && iface.hasOriginalProperty(ImagesData::DeflectionComponent2, s));

		// TODO: for now this needs to be in arcsec, make this configurable
		const float *pOrigAx = iface.getOriginalProperties(ImagesData::DeflectionComponent1, s);
		const float *pOrigAy = iface.getOriginalProperties(ImagesData::DeflectionComponent2, s);
		const Vector2Df *pAlpha = iface.getAlphas(s);

		for (int i = 0 ; i < numPoints ; i++)
		{
			Vector2Df diff = pAlpha[i];
			diff *= unitScale;
			diff -= Vector2Df(pOrigAx[i], pOrigAy[i]);
			deflAngleDiffSum2 += diff.getLengthSquared();
		}
	}

	if (totalNumPts > 0)
		deflAngleDiffSum2 /= (float)totalNumPts;

	fitness = deflAngleDiffSum2;
	return true;
}

} // end namespace
