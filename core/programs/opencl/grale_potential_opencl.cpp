// Currently, the potential values that are calculated in the OpenCL code,
// are not really used anywhere. This is a helper program to get the values
// calculated by the OpenCL code, mainly for debugging
// This is largely copy-paste from grale_lensplane_opencl.cpp

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
	string getRenderType() const { return "POTENTAL"; }
	string getVersionInfo() const { return "OpenCL potential renderer"; }
private:
	bool_t renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
	                  double x0, double y0, double dX, double dY, int numX, int numY,
	                  vector<double> &renderPoints);
	bool_t renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                         const std::vector<double> &inputXY, std::vector<double> &renderPoints)
	{
		return "Not supported";
	}

	OpenCLKernel *m_pCl;
};

bool_t OpenCLRenderer::renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
                                  double x0, double y0, double dX, double dY, int hNumX, int hNumY,
                                  vector<double> &renderPoints)
{
	renderPoints.resize(hNumX*hNumY);

	int numIntParams, numFloatParams;
	double deflectionScale, potentialScale;

	if (!pLens->getCLParameterCounts(&numIntParams, &numFloatParams))
		return "Couldn't get parameter counts for OpenCL program: " + pLens->getErrorString();

	if (!pLens->getSuggestedScales(&deflectionScale, &potentialScale))
	{
		deflectionScale = ANGLE_ARCSEC;
		potentialScale = ANGLE_ARCSEC*ANGLE_ARCSEC;
	}

	vector<int> intParams(numIntParams+1); // +1 To avoid trying to upload 0 bytes
	vector<float> floatParams(numFloatParams+1);

	intParams[numIntParams] = -12345;
	floatParams[numFloatParams] = -12345.0f;

	if (!pLens->getCLParameters(deflectionScale, potentialScale, &(intParams[0]), &(floatParams[0])))
		return "Couldn't get parameters for OpenCL program: " + pLens->getErrorString();

	std::cerr << "Integer parameters (" << numIntParams << ")" << std::endl;
	for (auto x : intParams)
		std::cerr << "  " << x << std::endl;

	std::cerr << std::endl;
	std::cerr << "Float parameters (" << numFloatParams << ")" << std::endl;
	for (auto x : floatParams)
		std::cerr << "  " << x << std::endl;

	string program, failLog, subRoutine;

	program = pLens->getCLLensProgram(deflectionScale, potentialScale, subRoutine);

	program += "\n";
	program += "__kernel void renderPotential(float2 startCoord, float2 step, int numX, int numY,\n";
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
	program += "	pResults[i+j*numX] = r.potential;\n";
	program += "}\n";

	cerr << program << endl;

	if (!m_pCl->loadKernel(program, "renderPotential", failLog))
	{
		if (!failLog.empty())
			cerr << failLog << endl;

		return "Couldn't load OpenCL kernel: " + m_pCl->getErrorString();
	}

	int numtotal = hNumX*hNumY;

	std::vector<float> hostResults(numtotal);
	cl_int err = 0;
	cl_int err1 = 0;
	cl_int err2 = 0;
	cl_int err3 = 0;
	cl_context context = m_pCl->getContext();

	memset(&(hostResults[0]), 0, sizeof(float)*hostResults.size());

	cl_mem clIntValues = m_pCl->clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(int) * (numIntParams+1), &(intParams[0]), &err1);
	cl_mem clFloatValues = m_pCl->clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float) * (numFloatParams+1), &(floatParams[0]), &err2);
	cl_mem clResults = m_pCl->clCreateBuffer(context, CL_MEM_WRITE_ONLY|CL_MEM_COPY_HOST_PTR, sizeof(float) * numtotal, &(hostResults[0]), &err3);

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
		err |= m_pCl->clEnqueueReadBuffer(commandQueue, clResults, true, 0, sizeof(float)*numtotal, &(hostResults[0]), 0, NULL, NULL);
		err |= m_pCl->clFinish(commandQueue);

		for (int i = 0 ; i < numtotal ; i++)
		{
			assert(i < renderPoints.size());

			renderPoints[i] = potentialScale * (double)hostResults[i];
		}
	}

	m_pCl->clReleaseMemObject(clIntValues);
	m_pCl->clReleaseMemObject(clFloatValues);
	m_pCl->clReleaseMemObject(clResults);

	if (err != CL_SUCCESS)
		return "Error executing rendering kernel";

	return true;
}

int main(int argc, char *argv[])
{
	string library = OpenCLKernel::getLibraryName();
	OpenCLKernel clKernel;
	if (!clKernel.loadLibrary(library))
	{
		cerr << "Unable to load OpenCL library: " << clKernel.getErrorString() << endl;
		return -1;
	}

	if (!clKernel.init())
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
