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

	double thetaRandomess = 0.1*ANGLE_ARCSEC;

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

	// TODO: also test with uncert enabled
	bool_t r;
	if (!(r = clDef.init(thetas, {}, intParams, floatParams, changeableParamIdx, 
	                     prog, subRoutName, 0, 0, originParams, numOriginParams,
						 recalcThetaInfo
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

		for (auto t : allTracedThetas)
			cout << t.getX() << "," << t.getY() << endl;

		throw runtime_error("TODO");

		double maxDiff = 0;
		double minDiff = numeric_limits<double>::max();
		auto recordDiff = [&maxDiff, &minDiff](double d)
		{
			if (d < minDiff)
				minDiff = d;
			if (d > maxDiff)
				maxDiff = d;

			return d;
		};

		size_t numGenomes = 1;
		for (size_t i = 0 ; i < numGenomes ; i++)
		{
			unique_ptr<GravitationalLens> newLens = lens.createLensFromCLFloatParams(deflScale, potScale, floatParams.data());
			for (size_t pt = 0 ; pt < thetas.size() ; pt++)
			{
				Vector2Dd theta(thetas[pt].getX(), thetas[pt].getY());
				theta *= deflScale;

				Vector2Dd alpha;
				newLens->getAlphaVector(theta, &alpha);
				alpha /= deflScale;

				Vector2Df gpuAlpha = allAlphas[i*thetas.size() + pt];

				auto diff = [&recordDiff](Vector2Dd a, Vector2Df b)
				{
					Vector2Dd b2(b.getX(), b.getY());
					Vector2Dd diff = a;
					diff -= b2;
					return recordDiff(diff.getLength());
				};

				cerr << "Genome " << i << " point " << pt << endl;
				cerr << "Alpha: CPU:  " << alpha.getX() << "," << alpha.getY() << endl;
				cerr << "       GPU:  " << gpuAlpha.getX() << "," << gpuAlpha.getY() << endl;
				cerr << "       diff: " << diff(alpha, gpuAlpha) << endl;

				double axx, ayy, axy;
				newLens->getAlphaVectorDerivatives(theta, axx, ayy, axy);
				float gpuAxx = allAxx[i*thetas.size() + pt];
				float gpuAyy = allAyy[i*thetas.size() + pt];
				float gpuAxy = allAxy[i*thetas.size() + pt];

				auto diff3 = [&recordDiff](double x0, double y0, double z0, float x1, float y1, float z1)
				{
					x0 -= x1; y0 -= y1; z0 -= z1;
					return recordDiff(std::sqrt(x0*x0 + y0*y0 + z0*z0));
				};

				cerr << "Deriv: CPU:  " << axx << "," << ayy << "," << axy << endl;
				cerr << "       GPU:  " << gpuAxx << "," << gpuAyy << "," << gpuAxy << endl;
				cerr << "       diff: " << diff3(axx, ayy, axy, gpuAxx, gpuAyy, gpuAxy) << endl;

				// Fot the potential we need to look at differences
				Vector2Dd theta0(thetas[0].getX(), thetas[0].getY());
				theta0 *= deflScale;

				double pot0 = numeric_limits<double>::quiet_NaN();
				newLens->getProjectedPotential(1.0, 1.0, theta0, &pot0);
				pot0 /= potScale;

				double pot = numeric_limits<double>::quiet_NaN();
				newLens->getProjectedPotential(1.0, 1.0, theta, &pot);
				pot /= potScale;

				// Look at potential differences!
				pot -= pot0;

				float gpuPot = allPotentials[i*thetas.size() + pt];
				gpuPot -= allPotentials[i*thetas.size() + 0];

				auto diff1 = [&recordDiff](double p0, float p1)
				{
					return recordDiff(std::abs(p0 - (double)p1));
				};

				cerr << "Poten: CPU: " << pot << endl;
				cerr << "       GPU: " << gpuPot << endl;
				cerr << "       diff: " << diff1(pot, gpuPot) << endl;

				cerr << endl;
			}
		}
		cerr << "Min diff: " << minDiff << endl;
		cerr << "Max diff: " << maxDiff << endl;
	}

	return 0;
}

