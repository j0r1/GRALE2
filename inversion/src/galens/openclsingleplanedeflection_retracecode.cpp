#include "openclsingleplanedeflection.h"
#include "retraceparameters.h"

using namespace std;
using namespace errut;

namespace grale
{

string getMultiStepCode(const string &functionName, const size_t numIterations, const string &lensRoutineName)
{
	return R"XYZ(

float2 )XYZ" + functionName + R"XYZ((const float2 thetaStart, const float2 betaTarget, const float dfrac, float *pBestBetaDiffSize,
						__global const int *pIntParams, __global const float *pFloatParams)
{
	int numIterations = )XYZ" + to_string(numIterations) + R"XYZ(;
	float2 theta = thetaStart;
	LensQuantities r = )XYZ" + lensRoutineName + R"XYZ((theta, pIntParams, pFloatParams);
	float2 betaCur = theta - dfrac*(float2)(r.alphaX, r.alphaY);
	float2 betaDiff = betaTarget - betaCur;
	float bestBetaDiffSize = sqrt(betaDiff.x*betaDiff.x + betaDiff.y*betaDiff.y);
	float2 bestRetraceTheta = theta;

	while (numIterations > 0)
	{
		float mxx = 1.0-dfrac*r.axx;
		float myy = 1.0-dfrac*r.ayy;
		float mxy = -dfrac*r.axy;

		// Get inverse magnification matrix and multiply with betaDiff
		float denom = mxx*myy - mxy*mxy;
		float txDiff = ( myy*betaDiff.x - mxy*betaDiff.y)/denom;
		float tyDiff = (-mxy*betaDiff.x + mxx*betaDiff.y)/denom;
		float2 thetaDiff = (float2)(txDiff, tyDiff);

		while (numIterations > 0)
		{
			numIterations--;

			r = )XYZ" + lensRoutineName + R"XYZ((theta+thetaDiff, pIntParams, pFloatParams);
			betaCur = (theta+thetaDiff) - dfrac*(float2)(r.alphaX, r.alphaY);
			betaDiff = betaTarget - betaCur;

			float betaDiffSize = sqrt(betaDiff.x*betaDiff.x + betaDiff.y*betaDiff.y);
			if (betaDiffSize < bestBetaDiffSize)
			{
				theta = theta + thetaDiff;
				bestBetaDiffSize = betaDiffSize;
				bestRetraceTheta = theta;
				break;
			}
			else
				thetaDiff *= (float2)(0.5, 0.5);
		}
	}

	*pBestBetaDiffSize = bestBetaDiffSize;
	return bestRetraceTheta;
}
)XYZ";
}

bool_t getReprojectSubroutineCode(const string &lensRoutineName, const SingleStepNewtonTraceParams &retraceParams, string &subCode)
{
	subCode = R"XYZ(

float2 findRetraceTheta(const float2 thetaStart, const float2 betaTarget, const float dfrac, float *pBestBetaDiffSize,
						__global const int *pIntParams, __global const float *pFloatParams)
{
	float2 theta = thetaStart;
	LensQuantities r = )XYZ" + lensRoutineName + R"XYZ((theta, pIntParams, pFloatParams);
	float2 betaCur = theta - dfrac*(float2)(r.alphaX, r.alphaY);
	float2 betaDiff = betaTarget - betaCur;

	float mxx = 1.0-dfrac*r.axx;
	float myy = 1.0-dfrac*r.ayy;
	float mxy = -dfrac*r.axy;

	// Get inverse magnification matrix and multiply with betaDiff
	float denom = mxx*myy - mxy*mxy;
	float txDiff = ( myy*betaDiff.x - mxy*betaDiff.y)/denom;
	float tyDiff = (-mxy*betaDiff.x + mxx*betaDiff.y)/denom;
	float2 thetaDiff = (float2)(txDiff, tyDiff);

	theta += thetaDiff;

	// TODO: disable this? Set pBestBetaDiffSize to NaN?
	r = )XYZ" + lensRoutineName + R"XYZ((theta, pIntParams, pFloatParams);

	betaCur = theta - dfrac*(float2)(r.alphaX, r.alphaY);
	betaDiff = betaTarget - betaCur;

	*pBestBetaDiffSize = sqrt(betaDiff.x*betaDiff.x + betaDiff.y*betaDiff.y);
	return theta;
}
)XYZ";
	
	return true;
}

bool_t getReprojectSubroutineCode(const string &lensRoutineName, const MultiStepNewtonTraceParams &retraceParams, string &subCode)
{
	size_t numIterations = retraceParams.getNumberOfEvaluations();
	subCode = getMultiStepCode("findRetraceTheta", numIterations, lensRoutineName);
	return true;
}

bool_t getReprojectSubroutineCode(const string &lensRoutineName, const ExpandedMultiStepNewtonTraceParams &retraceParams, string &subCode)
{
	ExpandedMultiStepNewtonTraceParams::Layout layout = retraceParams.getLayout();
	size_t numEvalsPerStartPosition = retraceParams.getNumberOfEvaluationsPerStartPosition();
	size_t maxGridSteps = retraceParams.getMaximumNumberOfGridSteps();
	double gridSpacing = retraceParams.getGridSpacing(); // Note: these are already converted to the angular units
	double acceptThreshold = retraceParams.getAcceptanceThreshold();

	subCode += getMultiStepCode("findRetraceTheta_singlepoint", numEvalsPerStartPosition, lensRoutineName);	
	subCode += R"XYZ(

float2 findRetraceTheta_level(int level, const float dxy,
		                const float2 thetaStart, const float2 betaTarget, const float dfrac, float *pBestBetaDiffSize,
						__global const int *pIntParams, __global const float *pFloatParams)
{
	if (level <= 1)
		return findRetraceTheta_singlepoint(thetaStart, betaTarget, dfrac, pBestBetaDiffSize, pIntParams, pFloatParams);

	float totalBestBetaDiff = INFINITY;
	float2 totalBestRetraceTheta = (float2)(INFINITY, INFINITY);

)XYZ";

	if (layout == ExpandedMultiStepNewtonTraceParams::FullGrid)
	{
		subCode += R"XYZ(
	for (int X = -(level-1) ; X <= (level-1) ; X++)
	{
		float2 theta = thetaStart + dxy * (float2)(X, -(level-1));
		float curBestBetaDiff = INFINITY;
		float2 curBestRetraceTheta = findRetraceTheta_singlepoint(theta, betaTarget, dfrac, &curBestBetaDiff, pIntParams, pFloatParams);
		if (curBestBetaDiff < totalBestBetaDiff)
		{
			totalBestBetaDiff = curBestBetaDiff;
			totalBestRetraceTheta = curBestRetraceTheta;
		}
	}

	for (int X = -(level-1) ; X <= (level-1) ; X++)
	{
		float2 theta = thetaStart + dxy * (float2)(X, +(level-1));
		float curBestBetaDiff = INFINITY;
		float2 curBestRetraceTheta = findRetraceTheta_singlepoint(theta, betaTarget, dfrac, &curBestBetaDiff, pIntParams, pFloatParams);
		if (curBestBetaDiff < totalBestBetaDiff)
		{
			totalBestBetaDiff = curBestBetaDiff;
			totalBestRetraceTheta = curBestRetraceTheta;
		}
	}

	for (int Y = -(level-1)+1 ; Y <= (level-1)-1 ; Y++)
	{
		float2 theta = thetaStart + dxy * (float2)(-(level-1), Y);
		float curBestBetaDiff = INFINITY;
		float2 curBestRetraceTheta = findRetraceTheta_singlepoint(theta, betaTarget, dfrac, &curBestBetaDiff, pIntParams, pFloatParams);
		if (curBestBetaDiff < totalBestBetaDiff)
		{
			totalBestBetaDiff = curBestBetaDiff;
			totalBestRetraceTheta = curBestRetraceTheta;
		}
	}

	for (int Y = -(level-1)+1 ; Y <= (level-1)-1 ; Y++)
	{
		float2 theta = thetaStart + dxy * (float2)(+(level-1), Y);
		float curBestBetaDiff = INFINITY;
		float2 curBestRetraceTheta = findRetraceTheta_singlepoint(theta, betaTarget, dfrac, &curBestBetaDiff, pIntParams, pFloatParams);
		if (curBestBetaDiff < totalBestBetaDiff)
		{
			totalBestBetaDiff = curBestBetaDiff;
			totalBestRetraceTheta = curBestRetraceTheta;
		}
	}
)XYZ";
	}
	else
	{
		if (layout == ExpandedMultiStepNewtonTraceParams::Square)
		{
			subCode += R"XYZ(
	const int dxs[] = { -1, -1, +1, +1 };
	const int dys[] = { -1, +1, -1, +1 };
	const int numPos = 4;
)XYZ";
		}
		else if (layout == ExpandedMultiStepNewtonTraceParams::Diamond)
		{
			subCode += R"XYZ(
	const int dxs[] = { 0, 0, -1, +1 };
	const int dys[] = { -1, +1, 0, 0 };
	const int numPos = 4;
)XYZ";
		}
		else if (layout == ExpandedMultiStepNewtonTraceParams::EightNeighbours)
		{
			subCode += R"XYZ(
	const int dxs[] = { -1, -1, +1, +1, 0, 0, -1, +1 };
	const int dys[] = { -1, +1, -1, +1, -1, +1, 0, 0 };
	const int numPos = 8;
)XYZ";
		}
		else
			return "Invalid layout type " + to_string((int)layout);

		subCode += R"XYZ(

	// Level is 1 for central point
	float diffScale = dxy*(float)(level-1);
	for (int i = 0 ; i < numPos ; i++)
	{
		float2 theta = thetaStart + diffScale * (float2)(dxs[i], dys[i]);
		float curBestBetaDiff = INFINITY;
		float2 curBestRetraceTheta = findRetraceTheta_singlepoint(theta, betaTarget, dfrac, &curBestBetaDiff, pIntParams, pFloatParams);
		if (curBestBetaDiff < totalBestBetaDiff)
		{
			totalBestBetaDiff = curBestBetaDiff;
			totalBestRetraceTheta = curBestRetraceTheta;
		}
	}
)XYZ";
	}

	subCode += R"XYZ(
	*pBestBetaDiffSize = totalBestBetaDiff;
	return totalBestRetraceTheta;
}

float2 findRetraceTheta(const float2 thetaStart, const float2 betaTarget, const float dfrac, float *pBestBetaDiffSize,
						__global const int *pIntParams, __global const float *pFloatParams)
{
	const float dxy = )XYZ" + to_string(gridSpacing) + R"XYZ(;
	const int maxLevel = )XYZ" + to_string(maxGridSteps) + R"XYZ(;
	const float betaDiffThreshold = )XYZ" + to_string(acceptThreshold) + R"XYZ(;

	float totalBestBetaDiff = INFINITY;
	float2 totalBestRetraceTheta = (float2)(INFINITY, INFINITY);

	for (int level = 1 ; level <= maxLevel ; level++)
	{
		float curBestBetaDiffSize = INFINITY;
		float2 curBestRetraceTheta = findRetraceTheta_level(level, dxy, thetaStart, betaTarget, dfrac, &curBestBetaDiffSize,
		                                                    pIntParams, pFloatParams);

		//printf("level %d, curBestBetaDiffSize %g\n", level, curBestBetaDiffSize);
		if (curBestBetaDiffSize <= betaDiffThreshold)
		{
			*pBestBetaDiffSize = curBestBetaDiffSize;
			return curBestRetraceTheta;
		}

		if (curBestBetaDiffSize < totalBestBetaDiff)
		{
			totalBestBetaDiff = curBestBetaDiffSize;
			totalBestRetraceTheta = curBestRetraceTheta;
		}
	}

	*pBestBetaDiffSize = totalBestBetaDiff;
	return totalBestRetraceTheta;
}
)XYZ";

	return true;
}

template <class T> bool_t callReprojectSub(const string &lensRoutineName, const TraceParameters &retraceParams, string &subCode, bool &handled)
{
	handled = false;

	if (dynamic_cast<const T *>(&retraceParams))
	{
		handled = true;
		return getReprojectSubroutineCode(lensRoutineName, static_cast<const T &>(retraceParams), subCode);
	}
	return true; // not handled, but no error
}

bool_t OpenCLSinglePlaneDeflection::getReprojectSubroutineCode(const string &lensRoutineName, const TraceParameters &retraceParams, string &subCode)
{
	if (dynamic_cast<const NoTraceParameters *>(&retraceParams))
		return "Internal error: getReprojectSubroutineCode should not be called without retrace parameters";

	bool handled = false;
	bool_t r;

	if (!(r = callReprojectSub<SingleStepNewtonTraceParams>(lensRoutineName, retraceParams, subCode, handled)))
		return r;
	if (handled)
		return true;

	if (!(r = callReprojectSub<MultiStepNewtonTraceParams>(lensRoutineName, retraceParams, subCode, handled)))
		return r;
	if (handled)
		return true;

	if (!(r = callReprojectSub<ExpandedMultiStepNewtonTraceParams>(lensRoutineName, retraceParams, subCode, handled)))
		return r;
	if (handled)
		return true;

	return "Internal error: unknown parameter type of trace parameters";



}

}
