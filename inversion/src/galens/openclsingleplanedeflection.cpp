#include "openclsingleplanedeflection.h"
#include "xoshiro128plus.h"
#include "opencl_xoshiro128plus.h"

using namespace std;
using namespace errut;

namespace grale
{

using namespace oclutils;

#define KERNELNUMBER_CALCULATE_DEFLECTION			0
#define KERNELNUMBER_FILLIN_CHANGED_PARAMS			1
#define KERNELNUMBER_CALCULATE_RANDOM_POINTOFFSETS	2
#define KERNELNUMBER_FETCH_ORIGIN_PARAMETERS		3

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
					   const std::vector<size_t> changeableParameterIndices,
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
					   int devIdx,
					   uint64_t initialUncertSeed,
					   const std::vector<std::pair<size_t, std::string>> &originParameters,
					   size_t numOriginParameters)
{
	if (m_init)
		return "Already initialized";
	
	auto cl = make_unique<OpenCLMultiKernel<4>>();

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
			return "No changeable parameters";

		vector<cl_int> clChangeableParams;
		for (auto x : changeableParameterIndices)
			clChangeableParams.push_back((cl_int)x);
		clChangeableParams.push_back(-67890); // sentinel and avoid empty vector

		if (!(r = m_clChangeableParamIndices.realloc(*cl, ctx, sizeof(cl_int)*clChangeableParams.size())) ||
			!(r = m_clChangeableParamIndices.enqueueWriteBuffer(*cl, queue, clChangeableParams, true)))
			return "Can't copy changeable parameter indices to GPU: " + cleanup(r.getErrorString());

		// Make sure enough is allocated for use in getChangeableParametersFromOriginParameters
		if (!(r = m_clChangedParamsBuffer.realloc(*cl, ctx, sizeof(cl_float)*changeableParameterIndices.size() )))
			return "Can't allocate initial changed parameters buffer: " + r.getErrorString();

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
		// TODO add code for xoshiro128plus functions
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
			return "Origin parameters specified, but size was set to 0";

		if (originParameters.size() != changeableParameterIndices.size())
			return "Origin parameters array length must be same as changeable parameters";

		vector<bool> originParameterUsed(numOriginParameters, false);
		vector<cl_int> originParameterIndices;
		string transformationCode;
		string caseCode;

		for (size_t destIdx = 0 ; destIdx < originParameters.size() ; destIdx++)
		{
			const auto & [srcIdx, code ] = originParameters[destIdx];
			originParameterIndices.push_back(srcIdx);
			if (srcIdx >= originParameterUsed.size()) 
				return "Invalid source index " + to_string(srcIdx) + " for origin parameter (length is " + to_string(numOriginParameters) + ")";
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
				return "Index " + to_string(srcIdx) + " of the " + to_string(numOriginParameters) + " origin parameters is not used";
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
}
)XYZ";
		cerr << "Compiling:\n" << kernelCode << endl;
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
			return "No origin parameters specified, but expected size is not zero";
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

	auto queue = m_cl->getCommandQueue();
	bool_t r;
	if (!(r = m_clOriginParams.enqueueWriteBuffer(*m_cl, queue, originParams.data(), true)))
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

	setKernArg(0, sizeof(cl_int), (void*)&clNumChangeable);
	setKernArg(1, sizeof(cl_int), (void*)&clNumOrigin);
	setKernArg(2, sizeof(cl_int), (void*)&clNumParamSets);
	setKernArg(3, sizeof(cl_mem), (void*)&m_clOriginParams.m_pMem);
	setKernArg(4, sizeof(cl_mem), (void*)&m_clOriginParamIndices.m_pMem);
	setKernArg(5, sizeof(cl_mem), (void*)&m_clChangedParamsBuffer.m_pMem);

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

bool_t OpenCLSinglePlaneDeflection::calculateDeflection(const std::vector<float> &parameters,  // should have numChangebleParams * numParamSets length
									  std::vector<Vector2Df> &allAlphas,
									  std::vector<float> &allAxx,
									  std::vector<float> &allAyy,
									  std::vector<float> &allAxy,
									  std::vector<float> &allPotentials)
{
	if (!m_init)
		return "Not initialized";

	size_t changeableSize = m_changeableParameterIndices.size();
	
	size_t numParamSets = 0;
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

	bool_t r;
	auto ctx = m_cl->getContext();
	auto queue = m_cl->getCommandQueue();

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

		setKernArg(0, sizeof(cl_int), (void*)&clNumChangeable);
		setKernArg(1, sizeof(cl_int), (void*)&clNumOrigin);
		setKernArg(2, sizeof(cl_int), (void*)&clNumParamSets);
		setKernArg(3, sizeof(cl_mem), (void*)&m_clOriginParams.m_pMem);
		setKernArg(4, sizeof(cl_mem), (void*)&m_clOriginParamIndices.m_pMem);
		setKernArg(5, sizeof(cl_mem), (void*)&m_clChangedParamsBuffer.m_pMem);

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
					   int devIdx,
					   uint64_t initialUncertSeed,
					   const std::vector<std::pair<size_t, std::string>> &originParameters,
					   size_t numOriginParameters
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
							 deflectionKernelCode, lensRoutineName, devIdx, initialUncertSeed,
							 originParameters, numOriginParameters);
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
		bool_t r = calculateDeflection(m_floatBuffer, m_allAlphas, m_allAxx, m_allAyy, m_allAxy, m_allPotentials);
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

bool OpenCLSinglePlaneDeflectionInstance::getResultsForGenome(const eatk::FloatVectorGenome &genome,
							 vector<Vector2Df> &alphas, vector<float> &axx,
							 vector<float> &ayy, vector<float> &axy,
							 vector<float> &potential)
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

	assert(sizeof(Vector2Df) == sizeof(float)*2);
	memcpy(alphas.data(), m_allAlphas.data()+pointOffsetStart, m_numPoints*2*sizeof(float));
	memcpy(axx.data(), m_allAxx.data()+pointOffsetStart, m_numPoints*sizeof(float));
	memcpy(ayy.data(), m_allAyy.data()+pointOffsetStart, m_numPoints*sizeof(float));
	memcpy(axy.data(), m_allAxy.data()+pointOffsetStart, m_numPoints*sizeof(float));
	memcpy(potential.data(), m_allPotentials.data()+pointOffsetStart, m_numPoints*sizeof(float));

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
