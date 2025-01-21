#include "openclsingleplanedeflection.h"
#include "constants.h"
#include "plummerlens.h"
#include "randomnumbergenerator.h"
#include "lensplane.h"
#include "imageplane.h"
#include <iostream>
#include <stdexcept>
#include <limits>
#include <serut/vectorserializer.h>

using namespace grale;
using namespace std;
using namespace errut;

class MyLensPlane : public LensPlane
{
protected:
	void setFeedbackStatus(const std::string &msg) override
	{
		cout << "STATUS: " << msg << endl;
	}
	void setFeedbackPercentage(int pct) override
	{
//		cout << "PCT: " << pct << endl;
	}
};

int main(int argc, char *argv[])
{
	if (argc != 2)
		throw runtime_error("Specify a lens file");

	string errStr;
	unique_ptr<GravitationalLens> pLens;
	if (!GravitationalLens::load(argv[1], pLens, errStr))
		throw runtime_error(errStr);

	GravitationalLens &lens = *pLens;

	double deflScale = numeric_limits<double>::quiet_NaN();
	double potScale = numeric_limits<double>::quiet_NaN();

	if (!lens.getSuggestedScales(&deflScale, &potScale))
		throw runtime_error("Can't get scales: " + lens.getErrorString());

	cerr << "Deflection scale = " << deflScale/ANGLE_ARCSEC << " arcsec" << endl;
	cerr << "Potential scale = " << potScale/(ANGLE_ARCSEC*ANGLE_ARCSEC) << " arcsec^2" << endl;

	RandomNumberGenerator rng;

	MyLensPlane lp;
	serut::VectorSerializer vSer;
	lens.write(vSer);
	if (!lp.init(vSer, { -100*ANGLE_ARCSEC, -100*ANGLE_ARCSEC }, { 100*ANGLE_ARCSEC, 100*ANGLE_ARCSEC }, 512, 512))
		throw runtime_error("Can't init lens plane: " + lp.getErrorString());

	size_t numBetas = 10;
	size_t srcNum = 0;
	double betaRange = 20*ANGLE_ARCSEC;
	vector<Vector2Df> thetas;
	vector<pair<int, float>> recalcThetaInfo;

	double thetaRandomess = 2.1*ANGLE_ARCSEC;
	//double thetaRandomess = 0;

	while (numBetas > 0)
	{
		double x = (rng.getRandomDouble()*2.0*betaRange - betaRange);
		double y = (rng.getRandomDouble()*2.0*betaRange - betaRange);
		double dfrac = rng.getRandomDouble()*0.25 + 0.65;

		ImagePlane ip;
		if (!ip.init(&lp, 1, dfrac))
			throw runtime_error(ip.getErrorString());

		Vector2Dd beta(x, y);
		vector<Vector2Dd> thetaPoints;

		if (!ip.traceBeta(beta, thetaPoints))
			throw runtime_error("Can't trace beta: " + ip.getErrorString());

		if (thetaPoints.size() > 1)
		{
			for (auto pt : thetaPoints)
			{
				double rndX = (rng.getRandomDouble()-0.5)*thetaRandomess;
				double rndY = (rng.getRandomDouble()-0.5)*thetaRandomess;
				
				recalcThetaInfo.push_back({ srcNum, dfrac });
				thetas.push_back({ (float)((pt.getX() + rndX)/deflScale), (float)((pt.getY() + rndY)/deflScale) });
			}

			numBetas--;
			srcNum++;
		}
	}

	cerr << "Using " << srcNum << " sources with total of " << thetas.size() << " images" << endl;

	int numIntParams = 0, numFloatParams = 0;
	lens.getCLParameterCounts(&numIntParams, &numFloatParams);
	
	vector<int> intParams(numIntParams);
	vector<float> floatParams(numFloatParams);
	if (!lens.getCLParameters(deflScale, potScale, intParams.data(), floatParams.data()))
		throw runtime_error("Can't get OpenCL parameters: " + lens.getErrorString());

	string subRoutName;
	string prog = lens.getCLLensProgram(deflScale, potScale, subRoutName);
	
	vector<size_t> changeableParamIdx = { 0 };
	vector<float> changedParameters = { floatParams[0] };

	OpenCLSinglePlaneDeflection clDef;
	vector<pair<size_t, string>> originParams;
	size_t numOriginParams = 0;
	size_t numRetraceSteps = 5;

	//vector<float> inKernelThetaUncert;
	//double inKernelThetaUncertSize = 2*ANGLE_ARCSEC;
	//for (auto x : thetas)
	//	inKernelThetaUncert.push_back( (float)(inKernelThetaUncertSize/deflScale) );

	// TODO: also test with uncert enabled
	bool_t r;
	if (!(r = clDef.init(thetas, {}, intParams, floatParams, changeableParamIdx, 
	                     prog, subRoutName, 0, 0, originParams, numOriginParams,
						 recalcThetaInfo, numRetraceSteps
						 )))
		throw runtime_error("Can't init OpenCL calculation code: " + r.getErrorString());

	vector<Vector2Df> allAlphas;
	vector<float> allAxx, allAyy, allAxy;
	vector<float> allPotentials;
	vector<Vector2Df> allTracedThetas;
	vector<float> allSourcePlaneDiffs;

	{
		if (!(r = clDef.calculateDeflectionAndRetrace(changedParameters,
											allAlphas, allAxx, allAyy, allAxy, allPotentials,
											allTracedThetas, allSourcePlaneDiffs)))
			throw runtime_error("Can't calculate deflections: " + r.getErrorString());

		if (allTracedThetas.size() == 0 || allSourcePlaneDiffs.size() == 0)
			throw runtime_error("No trace info was deteced");

		assert(thetas.size() == allTracedThetas.size() && allTracedThetas.size() == allSourcePlaneDiffs.size());

		for (auto i = 0 ; i < thetas.size() ; i++)
		{
			auto tOrig = thetas[i];
			auto tNew = allTracedThetas[i];

			cout << tOrig.getX() << " " << tOrig.getY() << " " << tNew.getX() << " " << tNew.getY() << " " << allSourcePlaneDiffs[i] << endl;
		}
	}

	return 0;
}

