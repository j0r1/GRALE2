#include "openclsingleplanedeflection.h"
#include "xoshiro128plus.h"
#include "opencl_xoshiro128plus.h"
#include <map>

using namespace std;
using namespace errut;

namespace grale
{

using namespace oclutils;

#define KERNELNUMBER_CALCULATE_DEFLECTION			0
#define KERNELNUMBER_FILLIN_CHANGED_PARAMS			1
#define KERNELNUMBER_CALCULATE_RANDOM_POINTOFFSETS	2
#define KERNELNUMBER_FETCH_ORIGIN_PARAMETERS		3
#define KERNELNUMBER_BACKPROJECT_POINTS				4
#define KERNELNUMBER_AVERAGE_BETAS					5
#define KERNELNUMBER_REPROJECT_BETAS				6
#define KERNELNUMBER_NEGLOGPRIORS					7

OpenCLSinglePlaneDeflection::OpenCLSinglePlaneDeflection()
{
}

OpenCLSinglePlaneDeflection::~OpenCLSinglePlaneDeflection()
{
}

bool_t OpenCLSinglePlaneDeflection::init(const std::vector<Vector2Df> &thetas, // already transformed into the correct units
					   const std::vector<float> &thetaUncert, // may be empty, otherwise in correct units and same length as thetas 
					   const std::vector<int> &templateIntParameters, // these cannot change
					   const std::vector<float> &templateFloatParameters, // only floating point params can change
					   const std::vector<size_t> &changeableParameterIndices,
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
					   const std::string &extraClPriorCode,
					   int devIdx,
					   uint64_t initialUncertSeed,
					   const std::vector<std::pair<size_t, std::string>> &originParameters,
					   size_t numOriginParameters,
					   const std::vector<std::pair<int, float>> &recalcThetaInfo,
					   size_t numRetraceIterations
					   )
{
	if (m_init)
		return "Already initialized";
	
	auto cl = make_unique<OpenCLMultiKernel<NumKernels>>();

	string libName = cl->getLibraryName();
	if (!cl->loadLibrary(libName))
		return "Can't load library '" + libName + "': " + cl->getErrorString();

	if (devIdx < 0)
	{
		devIdx = cl->getRotatedDeviceIndex();
		if (devIdx < 0) // something went wrong
			return "Can't get rotated device index: " + cl->getErrorString();
	}
	
	if (!cl->init(devIdx))
		return "Can't init device " + to_string(devIdx) + ": " + cl->getErrorString();

	auto cleanup = [this, &cl](string r)
	{
		m_clThetas.dealloc(*cl);
		m_clIntParams.dealloc(*cl);
		m_clFloatParams.dealloc(*cl);
		m_clAllResults.dealloc(*cl);
		m_clChangedParamsBuffer.dealloc(*cl);
		m_clChangeableParamIndices.dealloc(*cl);
		m_clThetaUncerts.dealloc(*cl);
		m_clThetaWithAdditions.dealloc(*cl);
		m_clRngStates.dealloc(*cl);
		m_clOriginParamIndices.dealloc(*cl);
		m_clOriginParams.dealloc(*cl);
		return r;
	};

	auto ctx = cl->getContext();
	auto queue = cl->getCommandQueue();

	// Copy the thetas to GPU	
	vector<cl_float> thetasFlat;
	for (auto t: thetas)
	{
		thetasFlat.push_back(t.getX());
		thetasFlat.push_back(t.getY());
	}

	if (thetasFlat.size() == 0)
		return cleanup("No theta positions");

	bool_t r;
	if (!(r = m_clThetas.realloc(*cl, ctx, sizeof(cl_float)*thetasFlat.size())) ||
		!(r = m_clThetas.enqueueWriteBuffer(*cl, queue, thetasFlat, true)))
		return "Can't copy thetas to GPU: " + cleanup(r.getErrorString());

	if (thetaUncert.size() > 0)
	{
		if (thetaUncert.size() != thetas.size())
			return cleanup("Must be same amount of uncertainties as positions");
		if (initialUncertSeed == 0)
			return cleanup("Initial seed for random position offsets must not be zero");
	
		vector<cl_float> clUncert;
		for (auto x : thetaUncert)
			clUncert.push_back(x);

		if (!(r = m_clThetaUncerts.realloc(*cl, ctx, sizeof(cl_float)*clUncert.size())) ||
			!(r = m_clThetaUncerts.enqueueWriteBuffer(*cl, queue, clUncert, true)))
			return "Can't copy theta uncertainties to GPU: " + cleanup(r.getErrorString());
		if (sizeof(uint32_t) != sizeof(cl_uint))
			return cleanup("Internal error: cl_uint doesn't appeat to be 32 bits?");

		vector<uint32_t> rngStates(clUncert.size()*4); // 4 uints for each point
		// Start first state from seed
		rngStates[0] = (uint32_t)((initialUncertSeed >> 0) & 0xffff);
		rngStates[1] = (uint32_t)((initialUncertSeed >> 16) & 0xffff);
		rngStates[2] = (uint32_t)((initialUncertSeed >> 32) & 0xffff);
		rngStates[3] = (uint32_t)((initialUncertSeed >> 48) & 0xffff);
		xoshiro128plus::jump(rngStates.data());

		// Set states for next points
		for (size_t i = 1 ; i < clUncert.size() ; i++)
			xoshiro128plus::jump(&rngStates[(i-1)*4], &rngStates[i*4]);

		//for (size_t i = 0 ; i < rngStates.size() ; i += 4)
		//	cerr << "rng state[" << (i/4) << "] = " << rngStates[i] << "," << rngStates[i+1] << "," << rngStates[i+2] << "," << rngStates[i+3] << endl;

		if (!(r = m_clRngStates.realloc(*cl, ctx, sizeof(uint32_t)*rngStates.size())) ||
			!(r = m_clRngStates.enqueueWriteBuffer(*cl, queue, rngStates, true)))
			return "Error uploading rng states: " + cleanup(r.getErrorString());

		vector<cl_float> thetaWithAdditions(clUncert.size()*2, numeric_limits<float>::quiet_NaN()); // initialize to NaN, *2 for two components
		if (!(r = m_clThetaWithAdditions.realloc(*cl, ctx, sizeof(cl_float)*thetaWithAdditions.size())) ||
			!(r = m_clThetaWithAdditions.enqueueWriteBuffer(*cl, queue, thetaWithAdditions, true)))
			return "Can't initialize theta additions to zero: " + cleanup(r.getErrorString());
	}
	else
	{
		if (initialUncertSeed != 0)
			return cleanup("Initial seed for random position offsets must be zero, since no randomization is requested");
	}

	vector<cl_int> clIntParams;
	for (auto x : templateIntParameters)
		clIntParams.push_back(x);
	clIntParams.push_back(-45678); // sentinel, and make sure that something is present
	
	vector<cl_float> clFloatParams;
	for (auto x : templateFloatParameters)
		clFloatParams.push_back(x);

	if (!(r = m_clIntParams.realloc(*cl, ctx, sizeof(cl_int)*clIntParams.size())) ||
		!(r = m_clIntParams.enqueueWriteBuffer(*cl, queue, clIntParams, true)))
		return "Can't copy integer parameters to GPU: " + cleanup(r.getErrorString());

	// Sanity checks
	if (changeableParameterIndices.size() > templateFloatParameters.size())
		return cleanup("Too many changeable indices");

	vector<bool> indexCheck(templateFloatParameters.size(), false);
	for (auto x : changeableParameterIndices)
	{
		if (x >= templateFloatParameters.size())
			return cleanup("Invalid index " + to_string(x) + " in changeable parameter indices");
		if (indexCheck[x])
			return cleanup("Trying to change index " + to_string(x) + " multiple times");
		indexCheck[x] = true;
	}

	// compile kernel(s)

	string deflectionKernel = deflectionKernelCode;
	deflectionKernel += R"XYZ(

__kernel void calculateDeflectionAngles(int numPoints, int numParamSets, int numFloatParams,
									__global const int *pIntParams,
									__global const float *pFloatParamsBase,
									__global const float *pThetas,
									__global float *pAllResults
									)
{
	const int i = get_global_id(0);
	if (i >= numPoints)
		return;
	const int paramSet = get_global_id(1);
	if (paramSet >= numParamSets)
		return;

	float2 theta;
	int thetaOffset = i*2;
	theta.x = pThetas[thetaOffset+0];
	theta.y = pThetas[thetaOffset+1];

	__global const float *pFloatParams = pFloatParamsBase + paramSet*numFloatParams;
	LensQuantities r = )XYZ" + lensRoutineName + R"XYZ((theta, pIntParams, pFloatParams);

	__global float *pResults = pAllResults + 6*numPoints*paramSet; // for each genome, 6 values
	int resultOffset = 6*i; // 6 resulting values per point

	pResults[resultOffset + 0] = r.alphaX;
	pResults[resultOffset + 1] = r.alphaY;
	pResults[resultOffset + 2] = r.axx;
	pResults[resultOffset + 3] = r.ayy;
	pResults[resultOffset + 4] = r.axy;
	pResults[resultOffset + 5] = r.potential;
}
)XYZ";

	// cerr << "Compiling OpenCL program:" << endl;
	// cerr << deflectionKernel << endl;

	string faillog;
	if (!cl->loadKernel(deflectionKernel, "calculateDeflectionAngles", faillog, KERNELNUMBER_CALCULATE_DEFLECTION))
	{
		cerr << faillog << endl;
		return cleanup(cl->getErrorString());
	}

	// Upload only the changed parameters (we'll need an extra kernel)
	{
		if (changeableParameterIndices.size() == 0)
			return cleanup("No changeable parameters");

		vector<cl_int> clChangeableParams;
		for (auto x : changeableParameterIndices)
			clChangeableParams.push_back((cl_int)x);
		clChangeableParams.push_back(-67890); // sentinel and avoid empty vector

		if (!(r = m_clChangeableParamIndices.realloc(*cl, ctx, sizeof(cl_int)*clChangeableParams.size())) ||
			!(r = m_clChangeableParamIndices.enqueueWriteBuffer(*cl, queue, clChangeableParams, true)))
			return "Can't copy changeable parameter indices to GPU: " + cleanup(r.getErrorString());

		// Make sure enough is allocated for use in getChangeableParametersFromOriginParameters
		if (!(r = m_clChangedParamsBuffer.realloc(*cl, ctx, sizeof(cl_float)*changeableParameterIndices.size() )))
			return "Can't allocate initial changed parameters buffer: " + cleanup(r.getErrorString());

		// And create kernel to incorporate changed parameters into full parameters

		string src = R"XYZ(

__kernel void fillInChangedParameters(__global const float *pChangedParamsBase,
									  __global float *pFloatParamsBase,
									  __global const int *pIndices,
									  int numFloatParams,
									  int numChangedParams,
									  int numParamSets)
{
	const int i = get_global_id(0);
	if (i >= numChangedParams)
		return;
	const int paramSet = get_global_id(1);
	if (paramSet >= numParamSets)
		return;

	__global const float *pChangedParams = pChangedParamsBase + paramSet*numChangedParams;
	__global float *pFloatParams = pFloatParamsBase + paramSet*numFloatParams;
	int destIdx = pIndices[i];
	pFloatParams[destIdx] = pChangedParams[i];
}

)XYZ";

		if (!cl->loadKernel(src, "fillInChangedParameters", faillog, KERNELNUMBER_FILLIN_CHANGED_PARAMS))
		{
			cerr << faillog << endl;
			return cleanup(cl->getErrorString());
		}
	}

	if (thetaUncert.size() > 0)
	{
		// kernel for randomization: based on thetas, states, sigmas, calculate thetas with additions
		string src = getOpenCLXoshiro128plusCode();
		src += R"XYZ(
__kernel void randomizeImagePlanePositions(int numPoints, 
									  __global const float *pThetas,
									  __global const float *pThetaUncerts,
									  __global float *pThetasWithAdditions,
									  __global uint *pRngStates)
{
	const int i = get_global_id(0);
	if (i >= numPoints)
		return;

	float2 theta = (float2)(pThetas[i*2], pThetas[i*2+1]);
	float sigma = pThetaUncerts[i];

	if (sigma == 0) // TODO: if it's zero, it stays zero, perhaps we shouldn't keep setting this?
	{
		pThetasWithAdditions[i*2] = theta.x;
		pThetasWithAdditions[i*2+1] = theta.y;
	}
	else
	{
		float2 randomGaussians = xoshiro128plus_next_gaussians(pRngStates + i*4);
		randomGaussians *= sigma;

		pThetasWithAdditions[i*2] = theta.x + randomGaussians.x;
		pThetasWithAdditions[i*2+1] = theta.y + randomGaussians.y;
	}
}
)XYZ";
		//cerr << "Compiling:\n" << src << endl;
		if (!cl->loadKernel(src, "randomizeImagePlanePositions", faillog, KERNELNUMBER_CALCULATE_RANDOM_POINTOFFSETS))
		{
			cerr << faillog << endl;
			return cleanup(cl->getErrorString());
		}
	}

	// See if we need to use 'origin parameters', so some floating point
	// parameters can be the same, or have some other relationship
	if (originParameters.size() > 0)
	{
		if (numOriginParameters == 0)
			return cleanup("Origin parameters specified, but size was set to 0");

		if (originParameters.size() != changeableParameterIndices.size())
			return cleanup("Origin parameters array length must be same as changeable parameters");

		vector<bool> originParameterUsed(numOriginParameters, false);
		vector<cl_int> originParameterIndices;
		string transformationCode;
		string caseCode;

		for (size_t destIdx = 0 ; destIdx < originParameters.size() ; destIdx++)
		{
			const auto & [srcIdx, code ] = originParameters[destIdx];
			originParameterIndices.push_back(srcIdx);
			if (srcIdx >= originParameterUsed.size()) 
				return cleanup("Invalid source index " + to_string(srcIdx) + " for origin parameter (length is " + to_string(numOriginParameters) + ")");
			originParameterUsed[srcIdx] = true;

			if (code.length() > 0) // code should be in curly braces, with return statement
			{
				string functionName = "originparameter_" + to_string(destIdx) + "_transformation";
				transformationCode += "\nfloat " + functionName + "(float x)\n";
				transformationCode += code;

				caseCode += R"XYZ(
	case )XYZ" + to_string(destIdx) + R"XYZ(:
		transformedValue = )XYZ" + functionName + R"XYZ((originValue);
		break;
)XYZ";
			}
		}

		for (size_t srcIdx = 0 ; srcIdx < originParameterUsed.size() ; srcIdx++)
		{
			if (!originParameterUsed[srcIdx])
				return cleanup("Index " + to_string(srcIdx) + " of the " + to_string(numOriginParameters) + " origin parameters is not used");
		}

		if (!(r = m_clOriginParamIndices.realloc(*cl, ctx, sizeof(cl_int)*originParameterIndices.size())) ||
			!(r = m_clOriginParamIndices.enqueueWriteBuffer(*cl, queue, originParameterIndices, true)))
			return "Can't copy origin parameter indices to GPU: " + cleanup(r.getErrorString());

		// we'll need to realloc this later, but we'll make sure that there's enough space for
		// one parameter set
		if (!(r = m_clOriginParams.realloc(*cl, ctx, sizeof(cl_float)*numOriginParameters)))
			return "Can't allocate buffer for the origin parameters: " + cleanup(r.getErrorString());

		// make kernel that transforms the origin parameters to the changeable parameters
		string kernelCode = transformationCode + R"XYZ(
__kernel void fetchOriginParameters(int numChangeableParams, int numOriginParams, int numParamSets,
							   __global const float *pAllOriginParams,
							   __global const int *pIndexMapping,
							   __global float *pAllChangeableParams)
							   //int debug)
{
	const int dstIdx = get_global_id(0);
	if (dstIdx >= numChangeableParams)
		return;
	const int paramSet = get_global_id(1);
	if (paramSet >= numParamSets)
		return;

	__global const float *pOriginParams = pAllOriginParams + paramSet*numOriginParams;
	__global float *pChangeableParams = pAllChangeableParams + paramSet*numChangeableParams;

	int srcIdx = pIndexMapping[dstIdx];
	float originValue = pOriginParams[srcIdx];
	float transformedValue = 0;

	switch(dstIdx)
	{
)XYZ";
		kernelCode += caseCode + R"XYZ(
	default:
		 transformedValue = originValue;
	}
	pChangeableParams[dstIdx] = transformedValue;
	//if (debug)
	//	printf("dstIdx = %d, srcIdx = %d, originValue = %g, transformedValue = %g\n", dstIdx, srcIdx, originValue, transformedValue);
}
)XYZ";
		//cerr << "Compiling:\n" << kernelCode << endl;
		string faillog;
		if (!cl->loadKernel(kernelCode, "fetchOriginParameters", faillog, KERNELNUMBER_FETCH_ORIGIN_PARAMETERS))
		{
			cerr << faillog << endl;
			return cleanup(cl->getErrorString());
		}
	}
	else
	{
		if (numOriginParameters != 0)
			return cleanup("No origin parameters specified, but expected size is not zero");
	}

	// Prior kernel
	if (extraClPriorCode.length() > 0)
	{
		string clPriorKernel = extraClPriorCode;
		clPriorKernel += R"XYZ(

__kernel void calculateTotalNegLogProbContributionsKernel(int numParamSets, int numFloatParams,
									__global const float *pFloatParamsBase,
									__global float *pResults)
{
	const int paramSet = get_global_id(0); // Only parallellize over parameter sets (genomes)
	if (paramSet >= numParamSets)
		return;

	__global const float *pFloatParams = pFloatParamsBase + paramSet*numFloatParams;
	pResults[paramSet] = calculateTotalNegLogProbContributions(pFloatParams);
}
)XYZ";
		string faillog;
		if (!cl->loadKernel(clPriorKernel, "calculateTotalNegLogProbContributionsKernel", faillog, KERNELNUMBER_NEGLOGPRIORS))
		{
			cerr << faillog << endl;
			cerr << "// Entire kernel code is:" << endl;
			cerr << clPriorKernel << endl;
			return cleanup("Error 'neglogprob' kernel, more info on stderr (set inverters.debugDirectStderr to True)" + cl->getErrorString());
		}
	} 

	if (recalcThetaInfo.size() > 0)
	{
		if (recalcThetaInfo.size() != thetas.size())
			return cleanup("Recalculate thetas vector should be of equal length as the thetas vector (" + to_string(recalcThetaInfo.size())
					       + " != " + to_string(thetas.size()) + ")");
		if (!(r = initRecalc(thetas.size(), recalcThetaInfo, *cl, deflectionKernelCode, lensRoutineName, numRetraceIterations)))
			return cleanup(r.getErrorString());

		m_recalcThetas = true;
	}

	// Ok

	m_numOriginParams = numOriginParameters;
	m_floatParamsCopy = clFloatParams;
	m_numPoints = thetas.size();
	m_numFloatParams = templateFloatParameters.size();
	m_currentNumParamSets = 0;
	m_maxNumParamSets = 0;
	m_changeableParameterIndices = changeableParameterIndices;
	m_cl = move(cl);
	m_init = true;
	m_devIdx = devIdx;
	m_haveClPriors = (extraClPriorCode.length() > 0)?true:false;

	return true;
}

errut::bool_t OpenCLSinglePlaneDeflection::initRecalc(size_t numTotalPoints, const std::vector<std::pair<int, float>> &recalcThetaInfo,
													  OpenCLMultiKernel<NumKernels> &cl, const std::string &deflectionKernelCode,
													  const std::string &lensRoutineName, size_t numRetraceIterations)
{
	// Build maps for distance fractions and points for sources

	map<size_t, float> sourceDfrac;
	map<size_t, vector<size_t>> sourceThetas;

	for (size_t tIdx = 0 ; tIdx < recalcThetaInfo.size() ; tIdx++)
	{
		auto [ srcIdx, dfrac ] = recalcThetaInfo[tIdx];
		if (srcIdx < 0) // signals not used
			continue;

		size_t src = (size_t)srcIdx;
		auto itFrac = sourceDfrac.find(src);
		if (itFrac == sourceDfrac.end())
			sourceDfrac[src] = dfrac;
		else
		{
			if (itFrac->second != dfrac)
				return "Incompatible distance fraction for source " + to_string(src) + ", point idx " + to_string(tIdx) + ": expecting "
					   + to_string(itFrac->second) + " but got " + to_string(dfrac);
		}

		auto itPoints = sourceThetas.find(src);
		if (itPoints == sourceThetas.end())
			sourceThetas[src] = { tIdx };
		else
		{
			auto &positions = itPoints->second;
			assert(positions.size() >= 1);
			assert(tIdx >= 1);
			assert(tIdx < numTotalPoints);

			if (positions.back()+1 != tIdx)
				return "Expecting points from same source to be consecutive: for source " + to_string(src) + " index " + to_string(tIdx)
					   + " does not follow " + to_string(positions.back());
			positions.push_back(tIdx);
		}
	}

	// Build mappings, in general these arrays will be shorter than the total number of points
	// The backproject kernel and reproject kernel will have this size in one dimension

	vector<cl_int> thetaIndices;
	vector<cl_float> imgDfracs;

	// Where the thetaIndices start for each source, and how many points it has for that source
	// This is needed to calculate the average of the backprojected points
	vector<cl_int> sourceStart;
	vector<cl_int> sourceNumImages;

	for (const auto& [ src, thetas ] : sourceThetas)
	{
		assert(sourceDfrac.find(src) != sourceDfrac.end());
		float dfrac = sourceDfrac[src];

		sourceStart.push_back((cl_int)thetaIndices.size()); // this is the start pos at that instant
		sourceNumImages.push_back((cl_int)thetas.size());

		for (auto t : thetas)
		{
			thetaIndices.push_back((cl_int)t);
			imgDfracs.push_back((cl_float)dfrac);
		}
	}

	assert(sourceStart.size() == sourceNumImages.size());
	m_clNumSources = (cl_int)sourceStart.size();

	assert(thetaIndices.size() == imgDfracs.size());
	assert(thetaIndices.size() <= numTotalPoints);
	m_clNumBpImages = (cl_int)thetaIndices.size();

	auto cleanup = [this, &cl](const string &r)
	{
		m_clBpDistFracs.dealloc(cl);
		m_clBpThetaIndices.dealloc(cl);
		m_clAllBetas.dealloc(cl);
		m_clSourceStarts.dealloc(cl);
		m_clSourceNumImages.dealloc(cl);
		m_clAllTracedThetas.dealloc(cl);
		m_clAllBetaDiffs.dealloc(cl);
		return r;
	};

	auto ctx = cl.getContext();
	auto queue = cl.getCommandQueue();
	bool_t r;

	if (!(r = m_clBpThetaIndices.realloc(cl, ctx, sizeof(cl_int)*thetaIndices.size())) ||
		!(r = m_clBpThetaIndices.enqueueWriteBuffer(cl, queue, thetaIndices, true)))
		return "Can't copy backproject theta indices to GPU: " + cleanup(r.getErrorString());

	if (!(r = m_clBpDistFracs.realloc(cl, ctx, sizeof(cl_float)*imgDfracs.size())) ||
		!(r = m_clBpDistFracs.enqueueWriteBuffer(cl, queue, imgDfracs, true)))
		return "Can't copy distance fractions to GPU: " + cleanup(r.getErrorString());

	// Backproject kernel
	{
		string backprojectKernel = deflectionKernelCode;
		backprojectKernel += R"XYZ(
__kernel void backprojectKernel(int numBpPoints, int numParamSets, int numFloatParams,
								__global const float *pFullThetas, // Can be either with or without randomization
								__global const int *pBpThetaIndices,
								__global const float *pDistFrac,
								__global const int *pIntParams,
								__global const float *pFloatParamsBase,
								__global float *pAllBetas
)
{
	const int ptIdx = get_global_id(0);
	if (ptIdx >= numBpPoints)
		return;
	const int paramSet = get_global_id(1);
	if (paramSet >= numParamSets)
		return;
	
	int thetaIdx = pBpThetaIndices[ptIdx] * 2; // *2 for two float components
	float2 theta = (float2)(pFullThetas[thetaIdx+0], pFullThetas[thetaIdx+1]);
	float dfrac = pDistFrac[ptIdx];

	__global const float *pFloatParams = pFloatParamsBase + paramSet*numFloatParams;
	LensQuantities r = )XYZ" + lensRoutineName + R"XYZ((theta, pIntParams, pFloatParams);

	__global float *pBetas = pAllBetas + 2*numBpPoints*paramSet;
	int resultOffset = 2*ptIdx;

	pBetas[resultOffset+0] = theta.x - dfrac * r.alphaX;
	pBetas[resultOffset+1] = theta.y - dfrac * r.alphaY;
}
)XYZ";
		string faillog;
		if (!cl.loadKernel(backprojectKernel, "backprojectKernel", faillog, KERNELNUMBER_BACKPROJECT_POINTS))
		{
			cerr << faillog << endl;
			return cleanup(cl.getErrorString());
		}
	}

	if (!(r = m_clSourceStarts.realloc(cl, ctx, sizeof(cl_int)*sourceStart.size())) ||
		!(r = m_clSourceStarts.enqueueWriteBuffer(cl, queue, sourceStart, true)))
		return "Can't copy source start indices to GPU: " + cleanup(r.getErrorString());

	if (!(r = m_clSourceNumImages.realloc(cl, ctx, sizeof(cl_int)*sourceNumImages.size())) ||
		!(r = m_clSourceNumImages.enqueueWriteBuffer(cl, queue, sourceNumImages, true)))
		return "Can't copy number of images for each source to GPU: " + cleanup(r.getErrorString());


	// Average source pos kernel
	{
		string averageBetaKernel = R"XYZ(
__kernel void averageSourcePositions(int numSources, int numBpPoints, int numParamSets,
									 __global const int *pSourceStarts,
									 __global const int *pSourceNumImages,
                                     __global float *pAllBetas
									 )
{
	const int srcIdx = get_global_id(0);
	if (srcIdx >= numSources)
		return;
	const int paramSet = get_global_id(1);
	if (paramSet >= numParamSets)
		return;

	__global float *pBetas = pAllBetas + 2*numBpPoints*paramSet;

	const int betaStartIdx = pSourceStarts[srcIdx] * 2;
	const int numPointsToAverage = pSourceNumImages[srcIdx];
	float2 srcAvg = (float2)(0.0, 0.0);

	for (int i = 0 ; i < numPointsToAverage ; i++)
		srcAvg += (float2)(pBetas[betaStartIdx + 2*i+0], pBetas[betaStartIdx + 2*i+1]);

	srcAvg /= (float)numPointsToAverage;

	// Store the average position for each image point
	for (int i = 0 ; i < numPointsToAverage ; i++)
	{
		pBetas[betaStartIdx + 2*i + 0] = srcAvg.x;
		pBetas[betaStartIdx + 2*i + 1] = srcAvg.y;
	}
}
)XYZ";
		string faillog;
		if (!cl.loadKernel(averageBetaKernel, "averageSourcePositions", faillog, KERNELNUMBER_AVERAGE_BETAS))
		{
			cerr << faillog << endl;
			cerr << averageBetaKernel << endl;
			return cleanup(cl.getErrorString());
		}
	}

	// Retrace kernel
	{
		string reprojectKernel = deflectionKernelCode;
		reprojectKernel += R"XYZ(
float2 retraceKernelStep(float2 theta, float dfrac, float2 betaTarget,
                         __global const int *pIntParams, __global const float *pFloatParams)
{
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
	return theta + (float2)(txDiff, tyDiff);
}

__kernel void retraceKernel(int numBpPoints, int numParamSets, int numFloatParams,
								__global const float *pFullThetas, // Can be either with or without randomization
								__global const int *pBpThetaIndices,
								__global const float *pDistFrac,
								__global const float *pAllAverageBetas,
								__global const int *pIntParams,
								__global const float *pFloatParamsBase,
								__global float *pAllTracedThetas,
								__global float *pAllSourcePlaneDists
)
{
	const int ptIdx = get_global_id(0);
	if (ptIdx >= numBpPoints)
		return;
	const int paramSet = get_global_id(1);
	if (paramSet >= numParamSets)
		return;

	int thetaIdx = pBpThetaIndices[ptIdx] * 2; // *2 for two float components
	float2 theta = (float2)(pFullThetas[thetaIdx+0], pFullThetas[thetaIdx+1]);
	float dfrac = pDistFrac[ptIdx];

	__global const float *pFloatParams = pFloatParamsBase + paramSet*numFloatParams;
	__global float *pAverageBetas = pAllAverageBetas + 2*numBpPoints*paramSet;

	const int resultOffset = 2*ptIdx;
	float2 betaTarget = (float2)(pAverageBetas[resultOffset + 0], pAverageBetas[resultOffset + 1]);

	// Do the refinement step for a number of iterations
	const int numIterations = )XYZ" + to_string(numRetraceIterations) + R"XYZ(;
	for (int i = 0 ; i < numIterations ; i++)
		theta = retraceKernelStep(theta, dfrac, betaTarget, pIntParams, pFloatParams);

	// And calculate the resulting difference in the source plane
	LensQuantities r = )XYZ" + lensRoutineName + R"XYZ((theta, pIntParams, pFloatParams);
	float2 betaCur = theta - dfrac*(float2)(r.alphaX, r.alphaY);
	float2 betaDiff = betaTarget - betaCur;
	float betaDiffSize = sqrt(betaDiff.x*betaDiff.x + betaDiff.y*betaDiff.y);

	__global float *pTracedThetas = pAllTracedThetas + 2*numBpPoints*paramSet;
	__global float *pSourcePlaneDists = pAllSourcePlaneDists + numBpPoints*paramSet;

	pTracedThetas[resultOffset+0] = theta.x;
	pTracedThetas[resultOffset+1] = theta.y;
	pSourcePlaneDists[ptIdx] = betaDiffSize;
}
)XYZ";
		string faillog;
		if (!cl.loadKernel(reprojectKernel, "retraceKernel", faillog, KERNELNUMBER_REPROJECT_BETAS))
		{
			cerr << faillog << endl;
			cerr << reprojectKernel << endl;
			return cleanup(cl.getErrorString());
		}
	}

	if (sizeof(Vector2Df) != sizeof(float)*2)
		return "Sanity check failed: 2D float vector is not of expected size";

	size_t testSize = 10;
	vector<Vector2Df> testVec(testSize);
	uint8_t *pStart = (uint8_t*)&(testVec[0]);
	uint8_t *pEnd = (uint8_t*)&(testVec[testSize-1]);
	size_t numBytes = pEnd-pStart;

	if (numBytes != 2*sizeof(float)*(testSize-1))
		return "Sanity check failed: vector<Vector2Df> is not laid out in memory as expected";

	return true;
}

void OpenCLSinglePlaneDeflection::destroy()
{
	if (!m_init)
		return;

	m_clThetas.dealloc(*m_cl);
	m_clIntParams.dealloc(*m_cl);
	m_clFloatParams.dealloc(*m_cl);
	m_clAllResults.dealloc(*m_cl);
	m_clChangedParamsBuffer.dealloc(*m_cl);
	m_clChangeableParamIndices.dealloc(*m_cl);
	m_clThetaUncerts.dealloc(*m_cl);
	m_clThetaWithAdditions.dealloc(*m_cl);
	m_clRngStates.dealloc(*m_cl);
	m_clOriginParamIndices.dealloc(*m_cl);
	m_clOriginParams.dealloc(*m_cl);

	m_clBpDistFracs.dealloc(*m_cl);
	m_clBpThetaIndices.dealloc(*m_cl);
	m_clAllBetas.dealloc(*m_cl);
	m_clSourceStarts.dealloc(*m_cl);
	m_clSourceNumImages.dealloc(*m_cl);
	m_clAllTracedThetas.dealloc(*m_cl);
	m_clAllBetaDiffs.dealloc(*m_cl);

	m_cl = nullptr;
	m_init = false;
}

bool_t OpenCLSinglePlaneDeflection::getChangeableParametersFromOriginParameters(const std::vector<float> &originParams,
															  std::vector<float> &changeableParams)
{
	if (!m_init)
		return "Not initialized";
	if (m_numOriginParams == 0)
		return "No origin parameters were specified at initialization";

	if (originParams.size() != m_numOriginParams)
		return "Incorrect number of origin parameters";

	changeableParams.resize(m_changeableParameterIndices.size());

	//cerr << "originParams = " << endl;
	//for (auto p : originParams)
	//	cerr << "   "  << p << endl;

	auto queue = m_cl->getCommandQueue();
	bool_t r;
	if (!(r = m_clOriginParams.enqueueWriteBuffer(*m_cl, queue, originParams, true)))
		return "Can't copy origin parameters to the GPU: " + r.getErrorString();

	auto kernel = m_cl->getKernel(KERNELNUMBER_FETCH_ORIGIN_PARAMETERS);
	cl_int err = 0;

	auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
	{
		err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
	};

	cl_int clNumChangeable = changeableParams.size();
	cl_int clNumOrigin = m_numOriginParams;
	cl_int clNumParamSets = 1;
	//cl_int debug = 1;

	setKernArg(0, sizeof(cl_int), (void*)&clNumChangeable);
	setKernArg(1, sizeof(cl_int), (void*)&clNumOrigin);
	setKernArg(2, sizeof(cl_int), (void*)&clNumParamSets);
	setKernArg(3, sizeof(cl_mem), (void*)&m_clOriginParams.m_pMem);
	setKernArg(4, sizeof(cl_mem), (void*)&m_clOriginParamIndices.m_pMem);
	setKernArg(5, sizeof(cl_mem), (void*)&m_clChangedParamsBuffer.m_pMem);
	//setKernArg(6, sizeof(cl_int), (void*)&debug);

	assert(m_clOriginParams.m_size >= sizeof(float)*clNumOrigin);
	assert(m_clChangedParamsBuffer.m_size >= sizeof(float)*clNumChangeable);

	if (err != CL_SUCCESS)
		return "Error setting kernel arguments for fetch origin parameters kernel";

	size_t workSize[2] = { changeableParams.size(), 1 };
	err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
	if (err != CL_SUCCESS)
		return "Error enqueing fetch origin parameter kernel: " + to_string(err);

	if (!(r = m_clChangedParamsBuffer.enqueueReadBuffer(*m_cl, queue, changeableParams, nullptr, nullptr, true)))
		return "Error reading resulting changable parameters: " + r.getErrorString();

	return true;
}

bool_t OpenCLSinglePlaneDeflection::getNumParamSets(const std::vector<float> &parameters, size_t &numParamSets)
{
	size_t changeableSize = m_changeableParameterIndices.size();

	numParamSets = 0;
	if (m_numOriginParams == 0) // Use changed parameters directly
	{
		if (parameters.size()%changeableSize != 0)
		   return "Invalid number of changeable parameters";
	
		numParamSets = parameters.size()/changeableSize;
	}
	else // The actual parameters are derived from some 'origin' parameters, possibly with a transformation
	{
		if (parameters.size()%m_numOriginParams != 0)
			return "Invalid number of origin parameters";

		numParamSets = parameters.size()/m_numOriginParams;
	}

	if (numParamSets == 0)
		return "No parameters specified";
	return true;
}

bool_t OpenCLSinglePlaneDeflection::calculateDeflection(const std::vector<float> &parameters,  // should have numChangebleParams * numParamSets length
									  std::vector<Vector2Df> &allAlphas,
									  std::vector<float> &allAxx,
									  std::vector<float> &allAyy,
									  std::vector<float> &allAxy,
									  std::vector<float> &allPotentials,
									  std::vector<float> &allClNegLogProbs)
{
	if (!m_init)
		return "Not initialized";

	bool_t r;
	size_t numParamSets = 0;
	if (!(r = getNumParamSets(parameters, numParamSets)))
		return r;

	auto ctx = m_cl->getContext();
	auto queue = m_cl->getCommandQueue();
	size_t changeableSize = m_changeableParameterIndices.size();

	// Need to do this here, in the following 'if' we'll be using the size
	m_allResultsBuffer.resize(numParamSets*m_numPoints*6); // 6 properties per point

	if (numParamSets > m_maxNumParamSets) // We need to prepare/upload some things
	{
		m_maxNumParamSets = numParamSets;

		vector<cl_float> clFloatParams;
		for (size_t i = 0 ; i < numParamSets ; i++)
			for (auto x : m_floatParamsCopy)
				clFloatParams.push_back(x);
		clFloatParams.push_back(-45678); // sentinel, and make sure that something is present

		if (!(r = m_clFloatParams.realloc(*m_cl, ctx, sizeof(cl_float)*clFloatParams.size())) ||
			!(r = m_clFloatParams.enqueueWriteBuffer(*m_cl, queue, clFloatParams, true)))
			return "Can't copy floating point parameters to GPU: " + r.getErrorString();

		swap(m_allFloatParams, clFloatParams);
		
		// allocate for outputs
		if (!(r = m_clAllResults.realloc(*m_cl, ctx, sizeof(cl_float)*m_allResultsBuffer.size())))
			return "Can't allocate results buffer: " + r.getErrorString();

		if (!(r = m_clChangedParamsBuffer.realloc(*m_cl, ctx, sizeof(cl_float)*(numParamSets*changeableSize + 1))))
			return "Can't allocate changed parameters buffer: " + r.getErrorString();

		if (m_numOriginParams > 0)
		{
			if (!(r = m_clOriginParams.realloc(*m_cl, ctx, sizeof(cl_float)*m_numOriginParams*numParamSets)))
				return "Can't allocate buffer for the origin parameters: " + r.getErrorString();
		}

		if (m_haveClPriors)
		{
			if (!(r = m_clPriorResults.realloc(*m_cl, ctx, sizeof(cl_float)*numParamSets)))
				return "Can't allocate buffer for the OpenCL prior results: " + r.getErrorString();
		}
	}

	m_currentNumParamSets = numParamSets;
	
	cl_int clNumParamSets = (cl_int)numParamSets;
	cl_int clNumFloatParams = (cl_int)m_numFloatParams;
	cl_int clNumChangeable = (cl_int)changeableSize;

	if (m_numOriginParams == 0)
	{
		// Upload the changed parameters to the GPU
		if (!(r = m_clChangedParamsBuffer.enqueueWriteBuffer(*m_cl, queue, parameters, true)))
			return "Can't copy changed parameters to GPU: " + r.getErrorString();
	}
	else // Get the real changed parameters from the specified 'origin' parameters
	{
		// Upload the origin parameters to the GPU an call a kernel to
		// transform these into the changed parameters

		if (!(r = m_clOriginParams.enqueueWriteBuffer(*m_cl, queue, parameters, true)))
			return "Can't copy origin parameters to the GPU: " + r.getErrorString();

		auto kernel = m_cl->getKernel(KERNELNUMBER_FETCH_ORIGIN_PARAMETERS);
		cl_int clNumOrigin = (cl_int)m_numOriginParams;
		cl_int err = 0;

		auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
		{
			err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
		};

		//cl_int debug = 0;

		setKernArg(0, sizeof(cl_int), (void*)&clNumChangeable);
		setKernArg(1, sizeof(cl_int), (void*)&clNumOrigin);
		setKernArg(2, sizeof(cl_int), (void*)&clNumParamSets);
		setKernArg(3, sizeof(cl_mem), (void*)&m_clOriginParams.m_pMem);
		setKernArg(4, sizeof(cl_mem), (void*)&m_clOriginParamIndices.m_pMem);
		setKernArg(5, sizeof(cl_mem), (void*)&m_clChangedParamsBuffer.m_pMem);
		//setKernArg(6, sizeof(cl_int), (void*)&debug);

		if (err != CL_SUCCESS)
			return "Error setting kernel arguments for fetch origin parameters kernel";

		size_t workSize[2] = { changeableSize, numParamSets };
		err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueing fetch origin parameter kernel: " + to_string(err);
	}

	// ... and fill in these changed parameters into the total parameters using a kernel
	{
		// trigger a kernel to fill in the changed parameters into the full ones
		auto kernel = m_cl->getKernel(KERNELNUMBER_FILLIN_CHANGED_PARAMS);
		cl_int err = 0;

		auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
		{
			err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
		};

		setKernArg(0, sizeof(cl_mem), (void*)&m_clChangedParamsBuffer.m_pMem);
		setKernArg(1, sizeof(cl_mem), (void*)&m_clFloatParams.m_pMem);
		setKernArg(2, sizeof(cl_mem), (void*)&m_clChangeableParamIndices.m_pMem);
		setKernArg(3, sizeof(cl_int), (void*)&clNumFloatParams);
		setKernArg(4, sizeof(cl_int), (void*)&clNumChangeable);
		setKernArg(5, sizeof(cl_int), (void*)&clNumParamSets);
		if (err != CL_SUCCESS)
			return "Error setting kernel arguments for parameter update kernel";

		size_t workSize[2] = { changeableSize, numParamSets };
		err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueing parameter fill in kernel: " + to_string(err);
	}

	// Prior code, if necessary
	if (m_haveClPriors)
	{
		auto kernel = m_cl->getKernel(KERNELNUMBER_NEGLOGPRIORS);
		cl_int err = 0;

		auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
		{
			err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
		};

		setKernArg(0, sizeof(cl_int), (void*)&clNumParamSets);
		setKernArg(1, sizeof(cl_int), (void*)&clNumFloatParams);
		setKernArg(2, sizeof(cl_mem), (void*)&m_clFloatParams.m_pMem);
		setKernArg(3, sizeof(cl_mem), (void*)&m_clPriorResults.m_pMem);

		if (err != CL_SUCCESS)
			return "Error setting kernel arguments for prior kernel";

		size_t workSize[1] = { numParamSets };
		err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, workSize, nullptr, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueing prior kernel: " + to_string(err);

		allClNegLogProbs.resize(numParamSets);
		if (!(r = m_clPriorResults.enqueueReadBuffer(*m_cl, queue, allClNegLogProbs, nullptr, nullptr, true)))
			return "Error reading prior results: " + r.getErrorString();
	}
	else
		allClNegLogProbs.clear();

	// Input position randomization is moved to other function

	// Finally perform the actual calculation of deflection angles

	auto kernel = m_cl->getKernel(KERNELNUMBER_CALCULATE_DEFLECTION);
	cl_int clNumPoints = (cl_int)m_numPoints;
	cl_int err = 0;

	auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
	{
		err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
	};

	// Depending on the case, either use the real thetas or the ones to which random numbers have been added
	void *pThetas = (m_clThetaWithAdditions.m_pMem)?m_clThetaWithAdditions.m_pMem:m_clThetas.m_pMem;

	setKernArg(0, sizeof(cl_int), (void*)&clNumPoints);
	setKernArg(1, sizeof(cl_int), (void*)&clNumParamSets);
	setKernArg(2, sizeof(cl_int), (void*)&clNumFloatParams);
	setKernArg(3, sizeof(cl_mem), (void*)&m_clIntParams.m_pMem);
	setKernArg(4, sizeof(cl_mem), (void*)&m_clFloatParams.m_pMem);
	setKernArg(5, sizeof(cl_mem), (void*)&pThetas);
	setKernArg(6, sizeof(cl_mem), (void*)&m_clAllResults.m_pMem);
	if (err != CL_SUCCESS)
		return "Error setting kernel arguments for deflection calculation kernel";

	size_t workSize[2] = { m_numPoints, numParamSets };
	err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
	if (err != CL_SUCCESS)
		return "Error enqueing calculation kernel: " + to_string(err);

	if (!(r = m_clAllResults.enqueueReadBuffer(*m_cl, queue, m_allResultsBuffer, nullptr, nullptr, true)))
		return "Error reading calculation results: " + r.getErrorString();

	// Copy the results to the arrays in the parameters
	size_t totalPoints = m_numPoints*numParamSets;
	allAlphas.resize(totalPoints);
	allAxx.resize(totalPoints);
	allAyy.resize(totalPoints);
	allAxy.resize(totalPoints);
	allPotentials.resize(totalPoints);

	for (size_t i = 0, j = 0 ; i < totalPoints ; i++, j += 6)
	{
		assert(j+5 < m_allResultsBuffer.size());
		allAlphas[i] = Vector2Df(m_allResultsBuffer[j+0], m_allResultsBuffer[j+1]);
		allAxx[i] = m_allResultsBuffer[j+2];
		allAyy[i] = m_allResultsBuffer[j+3];
		allAxy[i] = m_allResultsBuffer[j+4];
		allPotentials[i] = m_allResultsBuffer[j+5];
	}

	return true;
}


errut::bool_t OpenCLSinglePlaneDeflection::calculateDeflectionAndRetrace(const std::vector<float> &parameters,
								  std::vector<Vector2Df> &changedThetas,
								  std::vector<Vector2Df> &allAlphas,
								  std::vector<float> &allAxx,
								  std::vector<float> &allAyy,
								  std::vector<float> &allAxy,
								  std::vector<float> &allPotentials,
								  std::vector<float> &allClNegLogProbs,
								  std::vector<Vector2Df> &tracedThetas,
								  std::vector<float> &tracedBetaDiffs
								  )
{
	bool_t r;
	if (!(r = calculateDeflection(parameters, allAlphas, allAxx, allAyy, allAxy, allPotentials, allClNegLogProbs)))
		return r;

	auto ctx = m_cl->getContext();
	auto queue = m_cl->getCommandQueue();

	// Fetch the changed thetas
	if (m_clThetaWithAdditions.m_pMem)
	{
		assert(m_clThetaWithAdditions.m_size == sizeof(float)*2*m_numPoints);
		changedThetas.resize(m_numPoints);
		if (!(r = m_clThetaWithAdditions.enqueueReadBuffer(*m_cl, queue, changedThetas.data(), changedThetas.size()*sizeof(float)*2, nullptr, nullptr, true)))
			return "Can't copy randomized input positions: " + r.getErrorString();
	}
	else
	{
		changedThetas.clear();
	}

	if (!m_recalcThetas)
	{
		tracedThetas.clear();
		tracedBetaDiffs.clear();
		return true; // TODO: or generate error?
	}

	size_t numParamSets = 0;
	if (!(r = getNumParamSets(parameters, numParamSets)))
		return r;

	assert(numParamSets > 0);
	//cerr << "DEBUG: numParamSets = " << numParamSets << endl;

	// Allocate room to store (intermediate) betas
	if (!(r = m_clAllBetas.realloc(*m_cl, ctx, sizeof(cl_float)*2*m_clNumBpImages*numParamSets)))
		return "Can't allocate beta byffer: " + r.getErrorString();
	if (!(r = m_clAllTracedThetas.realloc(*m_cl, ctx, sizeof(cl_float)*2*m_clNumBpImages*numParamSets)))
		return "Can't alocate traced theta buffer on GPU: " + r.getErrorString();
	if (!(r = m_clAllBetaDiffs.realloc(*m_cl, ctx, sizeof(cl_float)*m_clNumBpImages*numParamSets)))
		return "Can't allocate source plane difference buffer on GPU: " + r.getErrorString();

	tracedThetas.resize(m_clNumBpImages*numParamSets);
	tracedBetaDiffs.resize(m_clNumBpImages*numParamSets);

	// Start first kernel
	{
		cl_int err = 0;
		auto kernel = m_cl->getKernel(KERNELNUMBER_BACKPROJECT_POINTS);
		auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
		{
			err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
		};

		cl_int clNumParamSets = (cl_int)numParamSets;
		cl_int clNumFloatParams = (cl_int)m_numFloatParams;
		// Depending on the case, either use the real thetas or the ones to which random numbers have been added
		void *pThetas = (m_clThetaWithAdditions.m_pMem)?m_clThetaWithAdditions.m_pMem:m_clThetas.m_pMem;

		setKernArg(0, sizeof(cl_int), (void*)&m_clNumBpImages);
		setKernArg(1, sizeof(cl_int), (void*)&clNumParamSets);
		setKernArg(2, sizeof(cl_int), (void*)&clNumFloatParams);
		setKernArg(3, sizeof(cl_mem), (void*)&pThetas);
		setKernArg(4, sizeof(cl_mem), (void*)&m_clBpThetaIndices.m_pMem);
		setKernArg(5, sizeof(cl_mem), (void*)&m_clBpDistFracs.m_pMem);
		setKernArg(6, sizeof(cl_mem), (void*)&m_clIntParams.m_pMem);
		setKernArg(7, sizeof(cl_mem), (void*)&m_clFloatParams.m_pMem);
		setKernArg(8, sizeof(cl_mem), (void*)&m_clAllBetas.m_pMem);

		if (err != CL_SUCCESS)
			return "Error setting kernel arguments for fetch origin parameters kernel";

		size_t workSize[2] = { (size_t)m_clNumBpImages, numParamSets };
		err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueing backproject parameter kernel: " + to_string(err);

		// For debugging! Disable this again!
		//if (!(r = m_clAllBetas.enqueueReadBuffer(*m_cl, queue, tracedThetas.data(), tracedThetas.size()*sizeof(float)*2, nullptr, nullptr, true)))
		//	return "Can't read betas: " + r.getErrorString();
	}

	// Average the betas
	{
		cl_int err = 0;
		auto kernel = m_cl->getKernel(KERNELNUMBER_AVERAGE_BETAS);
		auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
		{
			err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
		};

		cl_int clNumParamSets = (cl_int)numParamSets;

		setKernArg(0, sizeof(cl_int), (void*)&m_clNumSources);
		setKernArg(1, sizeof(cl_int), (void*)&m_clNumBpImages);
		setKernArg(2, sizeof(cl_int), (void*)&clNumParamSets);
		setKernArg(3, sizeof(cl_mem), (void*)&m_clSourceStarts);
		setKernArg(4, sizeof(cl_mem), (void*)&m_clSourceNumImages);
		setKernArg(5, sizeof(cl_mem), (void*)&m_clAllBetas.m_pMem);

		if (err != CL_SUCCESS)
			return "Error setting kernel arguments for fetch origin parameters kernel";

		size_t workSize[2] = { (size_t)m_clNumSources, numParamSets };
		err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueing averaging backprojected position kernel: " + to_string(err);

		// For debugging! Disable this again!
		//if (!(r = m_clAllBetas.enqueueReadBuffer(*m_cl, queue, tracedThetas.data(), tracedThetas.size()*sizeof(float)*2, nullptr, nullptr, true)))
		//	return "Can't read betas: " + r.getErrorString();
	}

	// Retrace the averaged beta positions
	{
		cl_int err = 0;
		auto kernel = m_cl->getKernel(KERNELNUMBER_REPROJECT_BETAS);
		auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
		{
			err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
		};

		cl_int clNumParamSets = (cl_int)numParamSets;
		cl_int clNumFloatParams = (cl_int)m_numFloatParams;
		// Depending on the case, either use the real thetas or the ones to which random numbers have been added
		void *pThetas = (m_clThetaWithAdditions.m_pMem)?m_clThetaWithAdditions.m_pMem:m_clThetas.m_pMem;

		setKernArg(0, sizeof(cl_int), (void*)&m_clNumBpImages);
		setKernArg(1, sizeof(cl_int), (void*)&clNumParamSets);
		setKernArg(2, sizeof(cl_int), (void*)&clNumFloatParams);
		setKernArg(3, sizeof(cl_mem), (void*)&pThetas);
		setKernArg(4, sizeof(cl_mem), (void*)&m_clBpThetaIndices.m_pMem);
		setKernArg(5, sizeof(cl_mem), (void*)&m_clBpDistFracs.m_pMem);
		setKernArg(6, sizeof(cl_mem), (void*)&m_clAllBetas.m_pMem);
		setKernArg(7, sizeof(cl_mem), (void*)&m_clIntParams.m_pMem);
		setKernArg(8, sizeof(cl_mem), (void*)&m_clFloatParams.m_pMem);
		setKernArg(9, sizeof(cl_mem), (void*)&m_clAllTracedThetas.m_pMem);
		setKernArg(10, sizeof(cl_mem), (void*)&m_clAllBetaDiffs.m_pMem);

		if (err != CL_SUCCESS)
			return "Error setting kernel arguments for fetch origin parameters kernel";

		size_t workSize[2] = { (size_t)m_clNumBpImages, numParamSets };
		err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 2, nullptr, workSize, nullptr, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueing retrace kernel: " + to_string(err);
		
		if (!(r = m_clAllTracedThetas.enqueueReadBuffer(*m_cl, queue, tracedThetas.data(), tracedThetas.size()*sizeof(float)*2, nullptr, nullptr, true)))
			return "Can't read betas: " + r.getErrorString();
		if (!(r = m_clAllBetaDiffs.enqueueReadBuffer(*m_cl, queue, tracedBetaDiffs, nullptr, nullptr, true)))
			return "Can't read source plane differences: " + r.getErrorString();
	}

	return true;
}

errut::bool_t OpenCLSinglePlaneDeflection::randomizeInputPositions()
{
	if (!m_init)
		return "Not initialized";

	if (m_clThetaWithAdditions.m_pMem) // Random additions are requested
	{
		auto ctx = m_cl->getContext();
		auto queue = m_cl->getCommandQueue();
		auto kernel = m_cl->getKernel(KERNELNUMBER_CALCULATE_RANDOM_POINTOFFSETS);
		cl_int clNumPoints = (cl_int)m_numPoints;
		cl_int err = 0;

		auto setKernArg = [this, &err, &kernel](cl_uint idx, size_t size, const void *ptr)
		{
			err |= m_cl->clSetKernelArg(kernel, idx, size, ptr);
		};

		setKernArg(0, sizeof(cl_int), (void*)&clNumPoints);
		setKernArg(1, sizeof(cl_mem), (void*)&m_clThetas.m_pMem);
		setKernArg(2, sizeof(cl_mem), (void*)&m_clThetaUncerts.m_pMem);
		setKernArg(3, sizeof(cl_mem), (void*)&m_clThetaWithAdditions.m_pMem);
		setKernArg(4, sizeof(cl_mem), (void*)&m_clRngStates);
		if (err != CL_SUCCESS)
			return "Error setting kernel arguments for randomization kernel";

		size_t workSize[1] = { m_numPoints };
		err = m_cl->clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, workSize, nullptr, 0, nullptr, nullptr);
		if (err != CL_SUCCESS)
			return "Error enqueing randomization kernel: " + to_string(err);
	}
	return true;
}


std::unique_ptr<OpenCLSinglePlaneDeflectionInstance> OpenCLSinglePlaneDeflectionInstance::s_instance;
std::mutex OpenCLSinglePlaneDeflectionInstance::s_instanceMutex;
std::set<uint64_t> OpenCLSinglePlaneDeflectionInstance::s_users;
bool OpenCLSinglePlaneDeflectionInstance::s_initTried = false;

errut::bool_t OpenCLSinglePlaneDeflectionInstance::initInstance(uint64_t userId,const std::vector<Vector2Df> &thetas, // already transformed into the correct units
					   const std::vector<float> &thetaUncert, // may be empty, otherwise in correct units and same length as thetas 
					   const std::vector<int> &templateIntParameters, // these cannot change
					   const std::vector<float> &templateFloatParameters, // only floating point params can change
					   const std::vector<size_t> changeableParameterIndices,
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
					   const std::string &extraClPriorCode,
					   int devIdx,
					   uint64_t initialUncertSeed,
					   const std::vector<std::pair<size_t, std::string>> &originParameters,
					   size_t numOriginParameters,
					   const std::vector<std::pair<int, float>> &recalcThetaInfo,
					   size_t numRetraceIterations
					   )

{
	lock_guard<mutex> guard(s_instanceMutex);

	if (s_users.find(userId) != s_users.end())
		return "An OpenCLSinglePlaneDeflectionInstance is already registered for this user";

	if (devIdx < 0)
		devIdx = -1;

	if (s_instance.get()) // Already initialized
	{
		if (devIdx != s_instance->getRequestedDeviceIndex())
			return "An OpenCLSinglePlaneDeflectionInstance already exists, but for a different device index (" + to_string(s_instance->getDeviceIndex()) + ") than requested (" + to_string(devIdx) + ")";

		cerr << "INFO: registered user in OpenCLSinglePlaneDeflectionInstance (2): " << userId << endl;
		s_users.insert(userId);
		return true;
	}

	if (s_initTried)
		return "GPU initialization failed before";
	s_initTried = true;

	unique_ptr<OpenCLSinglePlaneDeflectionInstance> oclCalc = make_unique<OpenCLSinglePlaneDeflectionInstance>();
	bool_t r = oclCalc->init(thetas, thetaUncert, templateIntParameters, templateFloatParameters, changeableParameterIndices,
							 deflectionKernelCode, lensRoutineName, extraClPriorCode, devIdx, initialUncertSeed,
							 originParameters, numOriginParameters, recalcThetaInfo, numRetraceIterations);
	oclCalc->m_requestedDevIdx = devIdx; // TODO: store this in a better way?

	if (!r)
		return "Can't init GPU: " + r.getErrorString();

	s_users.insert(userId);
	cerr << "INFO: registered user in OpenCLSinglePlaneDeflectionInstance: " << userId << endl;
	cerr << "INFO: requested GPU " << oclCalc->getRequestedDeviceIndex() << ", got " << oclCalc->getDeviceIndex() << endl;

	s_instance = move(oclCalc);
	return true;
}

void OpenCLSinglePlaneDeflectionInstance::releaseInstance(uint64_t userId)
{
	lock_guard<mutex> guard(s_instanceMutex);

	auto it = s_users.find(userId);
	if (it == s_users.end())
	{
		cerr << "WARNING: OpenCLSinglePlaneDeflectionInstance::releaseInstance: can't find user " << userId << endl;
		return;
	}
	s_users.erase(it);

	if (s_users.empty())
	{
		cerr << "INFO: no more users of OpenCLSinglePlaneDeflectionInstance, removing static instance" << endl;
		s_initTried = false;
		s_instance = nullptr;
	}
}

OpenCLSinglePlaneDeflectionInstance &OpenCLSinglePlaneDeflectionInstance::instance()
{
	auto pInst = s_instance.get();
	assert(pInst);
	return *pInst;
}

OpenCLSinglePlaneDeflectionInstance::OpenCLSinglePlaneDeflectionInstance()
{

}

OpenCLSinglePlaneDeflectionInstance::~OpenCLSinglePlaneDeflectionInstance()
{

}

// This is called on a new iteration, so that we know how much genomes are
// coming in before we're going to start the GPU calculation
void OpenCLSinglePlaneDeflectionInstance::setTotalGenomesToCalculate(size_t iteration, size_t num)
{
	lock_guard<mutex> guard(m_mutex);
	assert(num > 0);
	assert(m_totalGenomesToCalculate == 0 || m_totalGenomesToCalculate == num);

	m_totalGenomesToCalculate = num;

	// This function may be called by multiple threads, check when it's the first time
	if (iteration != m_prevIteration)
	{
		m_prevIteration = iteration;

		bool_t r = randomizeInputPositions();
		if (!r)
			cerr << "WARNING: Unexpected! Can't randomize input positions: " << r.getErrorString() << endl;
	}
}

bool_t OpenCLSinglePlaneDeflectionInstance::scheduleCalculation(const eatk::FloatVectorGenome &genome)
{
	lock_guard<mutex> guard(m_mutex);

	assert(!m_calculationDone);

	const auto *pGenome = &genome;
	assert(m_genomeOffsets.find(pGenome) == m_genomeOffsets.end()); // Make sure we don't have an entry yet

	const size_t numParamsToOptimize = (m_numOriginParams == 0)?m_changeableParameterIndices.size():m_numOriginParams;

	assert(genome.getValues().size() == numParamsToOptimize);
	size_t offset = m_floatBuffer.size()/numParamsToOptimize;
	for (auto x : genome.getValues())
		m_floatBuffer.push_back(x);

	m_genomeOffsets[pGenome] = offset; // Keep track where the results for this genome will be stored

	if (offset+1 == m_totalGenomesToCalculate)
	{
		//cerr << "Got " << m_totalGenomesToCalculate << " genomes, need to calculate!" << endl;
		bool_t r = calculateDeflectionAndRetrace(m_floatBuffer, m_adjustedThetas,
		                                         m_allAlphas, m_allAxx, m_allAyy, m_allAxy, m_allPotentials,
												 m_allClPriors,
												 m_allTracedThetas, m_allTracedBetaDiffs);
		if (!r)
			return r;
		//cerr << "Delections calculated successfully" << endl;
		m_calculationDone = true;

		// cerr << "GPU all alphas: " << endl;
		// for (auto a : m_allAlphas)
		//	 cerr << a.getX() << " " << a.getY() << endl;
	}
	return true;
}

bool OpenCLSinglePlaneDeflectionInstance::getAdjustedThetas(std::vector<Vector2Df> &adjustedThetas)
{
	lock_guard<mutex> guard(m_mutex);
	if (!m_calculationDone)
		return false;

	// We need to make a copy since multiple threads can request this information
	vector<Vector2Df> copy = m_adjustedThetas;
	adjustedThetas.swap(copy);

	return true;
}

bool OpenCLSinglePlaneDeflectionInstance::getResultsForGenome(const eatk::FloatVectorGenome &genome,
							 vector<Vector2Df> &alphas, vector<float> &axx,
							 vector<float> &ayy, vector<float> &axy,
							 vector<float> &potential,
							 vector<Vector2Df> &tracedThetas,
							 vector<float> &tracedBetaDiffs,
							 float &negLogPriorProb)
{
	lock_guard<mutex> guard(m_mutex);
	if (!m_calculationDone)
		return false;

	const auto *pGenome = &genome;
	auto it = m_genomeOffsets.find(pGenome);
	assert(it != m_genomeOffsets.end());

	int offset = it->second;
	int pointOffsetStart = offset*m_numPoints;
	int pointOffsetEnd = pointOffsetStart + m_numPoints;

	// cerr << "offset for genome " << (void*)pGenome << " is " << offset 
	//	  << ", numPoints = " << m_numPoints << endl;

	assert(pointOffsetStart < m_allAlphas.size());
	assert(pointOffsetEnd <= m_allAlphas.size());

	alphas.resize(m_numPoints);
	axx.resize(m_numPoints);
	ayy.resize(m_numPoints);
	axy.resize(m_numPoints);
	potential.resize(m_numPoints);

	tracedThetas.resize(m_clNumBpImages);
	tracedBetaDiffs.resize(m_clNumBpImages);

	assert(sizeof(Vector2Df) == sizeof(float)*2);
	memcpy(alphas.data(), m_allAlphas.data()+pointOffsetStart, m_numPoints*2*sizeof(float));
	memcpy(axx.data(), m_allAxx.data()+pointOffsetStart, m_numPoints*sizeof(float));
	memcpy(ayy.data(), m_allAyy.data()+pointOffsetStart, m_numPoints*sizeof(float));
	memcpy(axy.data(), m_allAxy.data()+pointOffsetStart, m_numPoints*sizeof(float));
	memcpy(potential.data(), m_allPotentials.data()+pointOffsetStart, m_numPoints*sizeof(float));

	if (m_clNumBpImages > 0)
	{
		int bpPointOffsetStart = offset*m_clNumBpImages;
		int bpPointOffsetEnd = bpPointOffsetStart + m_clNumBpImages;
		assert(bpPointOffsetStart < m_allTracedThetas.size());
		assert(bpPointOffsetStart <= m_allTracedThetas.size());

		memcpy(tracedThetas.data(), m_allTracedThetas.data() + bpPointOffsetStart, m_clNumBpImages*2*sizeof(float));
		memcpy(tracedBetaDiffs.data(), m_allTracedBetaDiffs.data() + bpPointOffsetStart, m_clNumBpImages*sizeof(float));
	}

	if (m_allClPriors.size() > 0)
	{
		assert(offset < m_allClPriors.size());
		negLogPriorProb = m_allClPriors[offset];
	}
	else
		negLogPriorProb = 0;

	m_genomeOffsets.erase(it);
	if (m_genomeOffsets.empty())
	{
		//cerr << "All results retrieved, ready to start again" << endl;
		m_calculationDone = false;
		m_totalGenomesToCalculate = 0;
		m_floatBuffer.resize(0);
	}

	return true;
}

bool_t OpenCLSinglePlaneDeflectionInstance::getChangeableParametersFromOriginParameters(const std::vector<float> &originParams,
												   std::vector<float> &changeableParams)
{
	lock_guard<mutex> guard(m_mutex);
	return OpenCLSinglePlaneDeflection::getChangeableParametersFromOriginParameters(originParams, changeableParams);
}

} // end namespace
