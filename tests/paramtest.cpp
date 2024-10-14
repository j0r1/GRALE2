#include "openclsingleplanedeflection.h"
#include "constants.h"
#include "plummerlens.h"
#include "randomnumbergenerator.h"
#include <iostream>
#include <stdexcept>
#include <limits>

using namespace grale;
using namespace std;
using namespace errut;

int main(void)
{
	PlummerLensParams params(1e13*MASS_SOLAR, 3*ANGLE_ARCSEC);
	PlummerLens lens;
	double Dd  = 1000*DIST_MPC;
	if (!lens.init(Dd, &params))
		throw runtime_error("Can't init lens: " + lens.getErrorString());

	double deflScale = numeric_limits<double>::quiet_NaN();
	double potScale = numeric_limits<double>::quiet_NaN();

	if (!lens.getSuggestedScales(&deflScale, &potScale))
		throw runtime_error("Can't get scales: " + lens.getErrorString());

	cerr << "Deflection scale = " << deflScale/ANGLE_ARCSEC << " arcsec" << endl;
	cerr << "Potential scale = " << potScale/(ANGLE_ARCSEC*ANGLE_ARCSEC) << " arcsec^2" << endl;

	RandomNumberGenerator rng;

	size_t numThetas = 20;
	double thetaRange = 100*ANGLE_ARCSEC;
	vector<Vector2Df> thetas(numThetas);
	for (auto &t : thetas)
	{
		double x = (rng.getRandomDouble()*2.0*thetaRange - thetaRange)/deflScale;
		double y = (rng.getRandomDouble()*2.0*thetaRange - thetaRange)/deflScale;
		t = Vector2Df((float)x, (float)y);
		cerr << " theta = " << x << ", " << y << endl;
	}

	int numIntParams = 0, numFloatParams = 0;
	lens.getCLParameterCounts(&numIntParams, &numFloatParams);
	
	vector<int> intParams(numIntParams);
	vector<float> floatParams(numFloatParams);
	if (!lens.getCLParameters(deflScale, potScale, intParams.data(), floatParams.data()))
		throw runtime_error("Can't get OpenCL parameters: " + lens.getErrorString());

	string subRoutName;
	string prog = lens.getCLLensProgram(deflScale, potScale, subRoutName);
	
	size_t numGenomes = 3;
	vector<size_t> changeableParamIdx = { 0, 1 };

	OpenCLSinglePlaneDeflection clDef;

	bool_t r;
	if (!(r = clDef.init(thetas, intParams, floatParams, changeableParamIdx, numGenomes, prog, subRoutName)))
		throw runtime_error("Can't init OpenCL calculation code: " + r.getErrorString());

	vector<Vector2Df> allAlphas;
	vector<float> allAxx, allAyy, allAxy;
	vector<float> allPotentials;

	vector<float> changedParams = {floatParams[0], floatParams[1],
	                               floatParams[0]/2, floatParams[1],
								   floatParams[0], floatParams[1]/2 };

	if (!(r = clDef.calculateDeflection(changedParams,
	                                    allAlphas, allAxx, allAyy, allAxy, allPotentials)))
		throw runtime_error("Can't calculate deflections: " + r.getErrorString());
	
	for (size_t i = 0 ; i < numGenomes ; i++)
	{
		vector<float> floatParamsMod = floatParams;
		for (size_t j = 0 ; j < changeableParamIdx.size() ; j++)
			floatParamsMod.at(changeableParamIdx[j]) = changedParams.at(i*changeableParamIdx.size() + j);
		
		unique_ptr<GravitationalLens> newLens = lens.createLensFromCLFloatParams(deflScale, potScale, floatParamsMod.data());
		for (size_t pt = 0 ; pt < thetas.size() ; pt++)
		{
			Vector2Dd theta(thetas[pt].getX(), thetas[pt].getY());
			theta *= deflScale;

			Vector2Dd alpha;
			newLens->getAlphaVector(theta, &alpha);
			alpha /= deflScale;

			Vector2Df gpuAlpha = allAlphas[i*thetas.size() + pt];

			auto diff = [](Vector2Dd a, Vector2Df b)
			{
				Vector2Dd b2(b.getX(), b.getY());
				Vector2Dd diff = a;
				diff -= b2;
				return diff.getLength();
			};

			cerr << "Genome " << i << " point " << pt << " CPU: " << alpha.getX() << "," << alpha.getY() 
			                                          << " GPU: " << gpuAlpha.getX() << "," << gpuAlpha.getY()
													  << " diff: " << diff(alpha, gpuAlpha) << endl;

			double axx, ayy, axy;
			newLens->getAlphaVectorDerivatives(theta, axx, ayy, axy);
			float gpuAxx = allAxx[i*thetas.size() + pt];
			float gpuAyy = allAyy[i*thetas.size() + pt];
			float gpuAxy = allAxy[i*thetas.size() + pt];

			auto diff3 = [](double x0, double y0, double z0, float x1, float y1, float z1)
			{
				x0 -= x1; y0 -= y1; z0 -= z1;
				return std::sqrt(x0*x0 + y0*y0 + z0*z0);
			};

			cerr << "Genome " << i << " point " << pt << " CPU: " << axx << "," << ayy << "," << axy
			                                          << " GPU: " << gpuAxx << "," << gpuAyy << "," << gpuAxy
													  << " diff: " << diff3(axx, ayy, axy, gpuAxx, gpuAyy, gpuAxy) << endl;

			cerr << endl;
		}
	}

	return 0;
}

