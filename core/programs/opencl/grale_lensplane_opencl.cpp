#include "communicator.h"
#include "gravitationallens.h"
#include "openclkernel.h"
#include "constants.h"
#include "utils.h"
#include <iostream>

using namespace std;
using namespace grale;
using namespace errut;

class OpenCLRenderer : public Communicator
{
public:
	OpenCLRenderer(OpenCLKernel *pClKernel) { m_pCl = pClKernel; }
	~OpenCLRenderer() { }
protected:
	string getRenderType() const { return "LENSPLANE"; }
	string getVersionInfo() const { return "OpenCL lensplane renderer"; }
private:
	bool_t renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
	                  double x0, double y0, double dX, double dY, int numX, int numY,
	                  vector<double> &renderPoints);
	bool_t renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                         const std::vector<double> &inputXY, std::vector<double> &renderPoints);

	OpenCLKernel *m_pCl;
};

bool_t OpenCLRenderer::renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
                                  double x0, double y0, double dX, double dY, int hNumX, int hNumY,
                                  vector<double> &renderPoints)
{
	renderPoints.resize(hNumX*hNumY*5);

	int numIntParams, numFloatParams;
	double deflectionScale, potentialScale;

	if (!pLens->getCLParameterCounts(&numIntParams, &numFloatParams))
		return "Couldn't get parameter counts for OpenCL program: " + pLens->getErrorString();

	if (!pLens->getSuggestedScales(&deflectionScale, &potentialScale))
	{
		deflectionScale = ANGLE_ARCSEC;
		potentialScale = ANGLE_ARCSEC*ANGLE_ARCSEC;
	}

	vector<int> intParams(numIntParams+1); // TODO: why did I do a +1 here?
	vector<float> floatParams(numFloatParams+1);

	if (!pLens->getCLParameters(deflectionScale, potentialScale, &(intParams[0]), &(floatParams[0])))
		return "Couldn't get parameters for OpenCL program: " + pLens->getErrorString();

	string program, failLog, subRoutine;

	program = pLens->getCLLensProgram(subRoutine);

	program += "\n";
	program += "__kernel void renderLensPlane(float2 startCoord, float2 step, int numX, int numY,\n";
	program += "                     __global const int *pIntParams,\n";
	program += "		     __global const float *pFloatParams,\n";
	program += "                     __global float *pResults)\n";
	program += "{\n";
	program += "	const int i = get_global_id(0);\n";
	program += "	const int j = get_global_id(1);\n";
	program += "\n";
	program += "	if (i >= numX || j >= numY)\n";
	program += "		return;\n";
	program += "\n";
	program += "	float2 pos = { startCoord.x + step.x*((float)i), startCoord.y + step.y*((float)j) } ;\n";
	program += "\n";
	program += "	LensQuantities r = " + subRoutine + "(pos, pIntParams, pFloatParams);\n";
	program += "\n";
	program += "	const int offset = (i+j*numX)*5;\n";
	program += "\n";
	program += "	pResults[offset+0] += r.alphaX;\n";
	program += "	pResults[offset+1] += r.alphaY;\n";
	program += "	pResults[offset+2] += r.axx;\n";
	program += "	pResults[offset+3] += r.ayy;\n";
	program += "	pResults[offset+4] += r.axy;\n";
	program += "}\n";

	cerr << program << endl;

	if (!m_pCl->loadKernel(program, "renderLensPlane", failLog))
	{
		if (!failLog.empty())
			cerr << failLog << endl;

		return "Couldn't load OpenCL kernel: " + m_pCl->getErrorString();

	}

	int numtotal = hNumX*hNumY;

	std::vector<float> hostResults(numtotal*5);
	cl_int err = 0;
	cl_int err1 = 0;
	cl_int err2 = 0;
	cl_int err3 = 0;
	cl_context context = m_pCl->getContext();

	memset(&(hostResults[0]), 0, sizeof(float)*hostResults.size());

	cl_mem clIntValues = m_pCl->clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(int) * (numIntParams+1), &(intParams[0]), &err1);
	cl_mem clFloatValues = m_pCl->clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float) * (numFloatParams+1), &(floatParams[0]), &err2);
	cl_mem clResults = m_pCl->clCreateBuffer(context, CL_MEM_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float) * numtotal * 5, &(hostResults[0]), &err3);

	err = err1|err2|err3;

	if (err != CL_SUCCESS)
	{
		if (clIntValues) m_pCl->clReleaseMemObject(clIntValues);
		if (clFloatValues) m_pCl->clReleaseMemObject(clFloatValues);
		if (clResults) m_pCl->clReleaseMemObject(clResults);

		return "Error allocating memory objects on GPU";
	}

	cl_int numX = hNumX;
	cl_int numY = hNumY;
	cl_int subStart = -1;
	cl_int subNum = -1;

	cl_float2 startPos = { (float)(x0/deflectionScale), (float)(y0/deflectionScale) };
	cl_float2 step = { (float)(dX/deflectionScale), (float)(dY/deflectionScale) };
	cl_kernel kernel = m_pCl->getKernel();

	err |= m_pCl->clSetKernelArg(kernel, 0, sizeof(cl_float2), (void *)&startPos);
	err |= m_pCl->clSetKernelArg(kernel, 1, sizeof(cl_float2), (void *)&step);
	err |= m_pCl->clSetKernelArg(kernel, 2, sizeof(cl_int), (void *)&numX);
	err |= m_pCl->clSetKernelArg(kernel, 3, sizeof(cl_int), (void *)&numY);
	err |= m_pCl->clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&clIntValues);
	err |= m_pCl->clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&clFloatValues);
	err |= m_pCl->clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&clResults);

	if (err != CL_SUCCESS)
	{
		m_pCl->clReleaseMemObject(clIntValues);
		m_pCl->clReleaseMemObject(clFloatValues);
		m_pCl->clReleaseMemObject(clResults);

		return "Error setting kernel arguments";
	}

	size_t workSize[2] = { (size_t)numX, (size_t)numY };
	cl_command_queue commandQueue = m_pCl->getCommandQueue();

	err |= m_pCl->clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL, workSize, NULL, 0, NULL, NULL);
	err |= m_pCl->clFinish(commandQueue);

	if (err == CL_SUCCESS)
	{
		err |= m_pCl->clEnqueueReadBuffer(commandQueue, clResults, true, 0, sizeof(float)*(numtotal*5), &(hostResults[0]), 0, NULL, NULL);
		err |= m_pCl->clFinish(commandQueue);

		int i5 = 0;

		for (int i = 0 ; i < numtotal ; i++, i5 += 5)
		{
			assert(i5+4 < renderPoints.size());

			renderPoints[i5+0] = deflectionScale * (double)hostResults[i5+0];
			renderPoints[i5+1] = deflectionScale * (double)hostResults[i5+1];
			renderPoints[i5+2] = hostResults[i5+2];
			renderPoints[i5+3] = hostResults[i5+3];
			renderPoints[i5+4] = hostResults[i5+4];
		}
	}

	m_pCl->clReleaseMemObject(clIntValues);
	m_pCl->clReleaseMemObject(clFloatValues);
	m_pCl->clReleaseMemObject(clResults);

	if (err != CL_SUCCESS)
		return "Error executing rendering kernel";

	return true;
}

bool_t OpenCLRenderer::renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
		                              const std::vector<double> &inputXY, std::vector<double> &renderPoints)
{
	int hNumXY = inputXY.size()/2;
	renderPoints.resize(hNumXY*5);

	int numIntParams, numFloatParams;
	double deflectionScale, potentialScale;

	if (!pLens->getCLParameterCounts(&numIntParams, &numFloatParams))
		return "Couldn't get parameter counts for OpenCL program: " + pLens->getErrorString();

	if (!pLens->getSuggestedScales(&deflectionScale, &potentialScale))
	{
		deflectionScale = ANGLE_ARCSEC;
		potentialScale = ANGLE_ARCSEC*ANGLE_ARCSEC;
	}

	vector<float> floatPos(inputXY.size());
	for (int i = 0 ; i < floatPos.size() ; i++)
		floatPos[i] = (float)(inputXY[i]/deflectionScale);

	vector<int> intParams(numIntParams+1); // TODO: why did I do a +1 here?
	vector<float> floatParams(numFloatParams+1);

	if (!pLens->getCLParameters(deflectionScale, potentialScale, &(intParams[0]), &(floatParams[0])))
		return "Couldn't get parameters for OpenCL program: " + pLens->getErrorString();

	string failLog, subRoutine;
	string program = pLens->getCLLensProgram(subRoutine);

	program += R"XYZ(
__kernel void renderLensPlane(int numXY, 
							__global const float *pInputPoints,
							__global const int *pIntParams,
							__global const float *pFloatParams,
							__global float *pResults)
{
	const int i = get_global_id(0);
	if (i >= numXY)
		return;

	float2 pos = { pInputPoints[i*2+0], pInputPoints[i*2+1] };

	LensQuantities r = )XYZ";

	program += subRoutine;
	program += R"XYZ((pos, pIntParams, pFloatParams);

	const int offset = i*5;

	pResults[offset+0] += r.alphaX;
	pResults[offset+1] += r.alphaY;
	pResults[offset+2] += r.axx;
	pResults[offset+3] += r.ayy;
	pResults[offset+4] += r.axy;
})XYZ";

	cerr << program << endl;

	if (!m_pCl->loadKernel(program, "renderLensPlane", failLog))
	{
		if (!failLog.empty())
			cerr << failLog << endl;

		return "Couldn't load OpenCL kernel: " + m_pCl->getErrorString();
	}

	std::vector<float> hostResults(hNumXY*5);
	cl_int err = 0;
	cl_int err0 = 0;
	cl_int err1 = 0;
	cl_int err2 = 0;
	cl_int err3 = 0;
	cl_context context = m_pCl->getContext();

	memset(&(hostResults[0]), 0, sizeof(float)*hostResults.size());

	cl_mem clFloatInputXY = m_pCl->clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float)*floatPos.size(), &floatPos[0], &err0);
	cl_mem clIntValues = m_pCl->clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(int) * (numIntParams+1), &(intParams[0]), &err1);
	cl_mem clFloatValues = m_pCl->clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float) * (numFloatParams+1), &(floatParams[0]), &err2);
	cl_mem clResults = m_pCl->clCreateBuffer(context, CL_MEM_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float) * hNumXY * 5, &(hostResults[0]), &err3);

	err = err0|err1|err2|err3;

	if (err != CL_SUCCESS)
	{
		if (clFloatInputXY) m_pCl->clReleaseMemObject(clFloatInputXY);
		if (clIntValues) m_pCl->clReleaseMemObject(clIntValues);
		if (clFloatValues) m_pCl->clReleaseMemObject(clFloatValues);
		if (clResults) m_pCl->clReleaseMemObject(clResults);

		return "Error allocating memory objects on GPU";
	}

	cl_int numXY = hNumXY;
	cl_int subStart = -1;
	cl_int subNum = -1;

	cl_kernel kernel = m_pCl->getKernel();

	err |= m_pCl->clSetKernelArg(kernel, 0, sizeof(cl_int), (void *)&numXY);
	err |= m_pCl->clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&clFloatInputXY);
	err |= m_pCl->clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&clIntValues);
	err |= m_pCl->clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&clFloatValues);
	err |= m_pCl->clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&clResults);

	if (err != CL_SUCCESS)
	{
		m_pCl->clReleaseMemObject(clFloatInputXY);
		m_pCl->clReleaseMemObject(clIntValues);
		m_pCl->clReleaseMemObject(clFloatValues);
		m_pCl->clReleaseMemObject(clResults);

		return "Error setting kernel arguments";
	}

	size_t workSize[1] = { (size_t)numXY };
	cl_command_queue commandQueue = m_pCl->getCommandQueue();

	err |= m_pCl->clEnqueueNDRangeKernel(commandQueue, kernel, 1, NULL, workSize, NULL, 0, NULL, NULL);
	err |= m_pCl->clFinish(commandQueue);

	if (err == CL_SUCCESS)
	{
		err |= m_pCl->clEnqueueReadBuffer(commandQueue, clResults, true, 0, sizeof(float)*(hNumXY*5), &(hostResults[0]), 0, NULL, NULL);
		err |= m_pCl->clFinish(commandQueue);

		int i5 = 0;

		for (int i = 0 ; i < hNumXY ; i++, i5 += 5)
		{
			assert(i5+4 < renderPoints.size());

			renderPoints[i5+0] = deflectionScale * (double)hostResults[i5+0];
			renderPoints[i5+1] = deflectionScale * (double)hostResults[i5+1];
			renderPoints[i5+2] = hostResults[i5+2];
			renderPoints[i5+3] = hostResults[i5+3];
			renderPoints[i5+4] = hostResults[i5+4];
		}
	}

	m_pCl->clReleaseMemObject(clFloatInputXY);
	m_pCl->clReleaseMemObject(clIntValues);
	m_pCl->clReleaseMemObject(clFloatValues);
	m_pCl->clReleaseMemObject(clResults);

	if (err != CL_SUCCESS)
		return "Error executing rendering kernel";

	return true;
}

int main(int argc, char *argv[])
{
#ifdef GRALE_LOADLIBRARY
	string defaultLibrary = "opencl.dll";
#else
	#ifdef __APPLE__
	string defaultLibrary = "/System/Library/Frameworks/OpenCL.framework/OpenCL";
#else
	string defaultLibrary = "libOpenCL.so";
#endif // __APPLE
#endif // GRALE_LOADLIBRARY

	string library = defaultLibrary;
	getenv("GRALE_OPENCLLIB", library); // Doesn't change the value if environment variable not set

	OpenCLKernel clKernel;
	if (!clKernel.init(library))
	{
		cerr << "Unable to initialize OpenCL: " << clKernel.getErrorString() << endl;
		return -1;
	}

	cerr << "OpenCL initialized" << endl;

	OpenCLRenderer renderer(&clKernel);
	bool_t r;

	if (!(r = renderer.render()))
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}
	return 0;
}
