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
	
	size_t numGenomes = 1;
	vector<size_t> changeableParamIdx = { 0 };

	OpenCLSinglePlaneDeflection clDef;

	bool_t r;
	if (!(r = clDef.init(thetas, intParams, floatParams, changeableParamIdx, numGenomes, prog )))
		throw runtime_error("Can't init OpenCL calculation code: " + r.getErrorString());
	return 0;
}
