#pragma once

#include "opencllibrary.h"
#include "pernodecounter.h"
#include "utils.h"
#include <vector>
#include <array>
#include <cassert>
#include <iostream>

template<int NumKernels>
class OpenCLMultiKernel : public OpenCLLibrary
{
public:
	OpenCLMultiKernel();
	~OpenCLMultiKernel();

	int getRotatedDeviceIndex();

	bool init(int devIdx = 0);
	bool loadKernel(const std::string &program, const std::string &kernelName, std::string &failLog, size_t kernelIdx = 0);
	void destroy();

	cl_context getContext()										{ return m_context; }
	cl_kernel getKernel(size_t kernelIdx = 0) 					{ assert(kernelIdx < NumKernels); return m_kernels[kernelIdx]; }
	cl_command_queue getCommandQueue()		 					{ return m_queue; }
private:
	void releaseAll();

	bool m_init;
	std::array<std::string, NumKernels> m_currentPrograms, m_currentKernels;
	
	cl_context m_context;
	cl_command_queue m_queue;
	std::array<cl_program, NumKernels> m_programs;
	std::array<cl_kernel, NumKernels> m_kernels;
	cl_device_id m_device;

	std::unique_ptr<grale::PerNodeCounter> m_perNodeCounter;
};

template<int NumKernels>
OpenCLMultiKernel<NumKernels>::OpenCLMultiKernel()
{
	m_init = false;
	m_context = 0;
	m_queue = 0;
	for (auto &p : m_programs)
		p = 0;
	for (auto &k : m_kernels)
		k = 0;
}

template<int NumKernels>
OpenCLMultiKernel<NumKernels>::~OpenCLMultiKernel()
{
	destroy();
}

template<int NumKernels>
bool OpenCLMultiKernel<NumKernels>::init(int devIdx)
{
	if (m_init)
	{
		setErrorString("OpenCL kernel loader is already initialized");
		return false;
	}

	for (auto &c : m_currentPrograms)
		c = "";
	for (auto &c : m_currentKernels)
		c = "";

	cl_platform_id platform;
	int numDevices, clDevType;

	if (!getPlatformAndDeviceCount(platform, numDevices, clDevType))
		return false;

	if (devIdx < 0 || devIdx >= numDevices)
	{
		setErrorString("Invalid device index (" + std::to_string(devIdx) + "), there are " + std::to_string(numDevices) + " detected");
		return false;
	}

	cl_int err;
	std::vector<cl_device_id> devices(numDevices);
	clGetDeviceIDs(platform, clDevType, numDevices, devices.data(), nullptr);

	std::string devTypeStr = "unknown";
	if (clDevType == CL_DEVICE_TYPE_CPU)
		devTypeStr = "CPU";
	else if (clDevType == CL_DEVICE_TYPE_GPU)
		devTypeStr = "GPU";
	
	std::cerr << "INFO: Found " << numDevices << " OpenCL devices of type " << devTypeStr << std::endl;
	std::cerr << "INFO: OpenCL Devices:" << std::endl;
	for (auto x : devices)
	{
		char name[1024] = { 0 };
		err = clGetDeviceInfo(x, CL_DEVICE_NAME, sizeof(name)-1, name, nullptr);
		if (err != CL_SUCCESS)
			std::cerr << "   (error " + std::to_string(err) + ")" <<std::endl;
		else
			std::cerr << "   " << x << ": " << name << std::endl;
	}

	m_device = devices[devIdx];
	m_context = clCreateContext(nullptr, 1, &m_device, nullptr, nullptr, &err);
	if (err != CL_SUCCESS)
	{
		setErrorString(std::string("Can't create OpenCL context: ") + getCLErrorString(err));
		releaseAll();
		return false;
	}

	m_queue = clCreateCommandQueue(m_context, m_device, 0, &err);
	if (err != CL_SUCCESS)
	{
		setErrorString(std::string("Can't create OpenCL command queue: ") + getCLErrorString(err));
		releaseAll();
		return false;
	}

#ifdef GRALE_DEBUG_OPENCLSENTINEL
	setSavedQueue(m_queue);
#endif // GRALE_DEBUG_OPENCLSENTINEL

	m_init = true;

	return true;
}

template<int NumKernels>
bool OpenCLMultiKernel<NumKernels>::loadKernel(const std::string &programString, const std::string &kernelName, std::string &failLog, size_t kernelIdx)
{
	assert(kernelIdx < NumKernels);

	if (!m_init)
	{
		setErrorString("OpenCL kernel loader was not initialized");
		return false;
	}

	if (programString.length() == 0 || kernelName.length() == 0)
	{
		setErrorString("No program or no kernel name specified");
		return false;
	}

	auto &program = m_programs[kernelIdx];
	auto &kernel = m_kernels[kernelIdx];
	auto &currentProgram = m_currentPrograms[kernelIdx];
	auto &currentKernel = m_currentKernels[kernelIdx];

	if (program == 0 || (program != 0 && programString != currentProgram))
	{
		// compile
		
		if (program)
			clReleaseProgram(program);
		if (kernel)
			clReleaseKernel(kernel);

		program = 0;
		kernel = 0;
		currentProgram = "";
		currentKernel = "";

		const char *pProgStr = programString.c_str();
		cl_int err;

		program = clCreateProgramWithSource(m_context, 1, &pProgStr, 0, &err);
		if (err != CL_SUCCESS)
		{
			setErrorString(std::string("Unable to create OpenCL program from specified source: ") + getCLErrorString(err));
			return false;
		}

		err = clBuildProgram(program, 1, &m_device, "-cl-std=CL2.0", nullptr, nullptr); // Needed CL2.0 support for printf for debugging
		//err = clBuildProgram(program, 1, &m_device, "", nullptr, nullptr);
		if (err != CL_SUCCESS)
		{
			size_t logLength;

			clGetProgramBuildInfo(program, m_device, CL_PROGRAM_BUILD_LOG, 0, 0, &logLength);

			std::vector<char> pLog(logLength+1);
			clGetProgramBuildInfo(program, m_device, CL_PROGRAM_BUILD_LOG, logLength, pLog.data(), 0);
			pLog[logLength] = 0;

			failLog = std::string(pLog.data());

			clReleaseProgram(program);
			program = 0;

			setErrorString(std::string("Unable to build OpenCL program: ") + getCLErrorString(err));

			return false;
		}

		currentProgram = programString;
	}

	if (kernel == 0 || (kernel != 0 && currentKernel != kernelName))
	{
		cl_int err;

		if (kernel)
			clReleaseKernel(kernel);

		currentKernel = "";

		kernel = clCreateKernel(program, kernelName.c_str(), &err);
		if (err != CL_SUCCESS)
		{
			setErrorString(std::string("Unable to get OpenCL kernel from program: ") + getCLErrorString(err));
			return false;
		}

		currentKernel = kernelName;
	}

	return true;
}

template<int NumKernels>
void OpenCLMultiKernel<NumKernels>::destroy()
{
	m_init = false;

	for (auto &c : m_currentKernels)
		c = "";
	for (auto &c : m_currentPrograms)
		c = "";

	releaseAll();
}

template<int NumKernels>
void OpenCLMultiKernel<NumKernels>::releaseAll()
{
	for (auto &k : m_kernels)
		if (k)
			clReleaseKernel(k);
	for (auto &p : m_programs)
		if (p)
			clReleaseProgram(p);

	if (m_queue)
		clReleaseCommandQueue(m_queue);
	if (m_context)
		clReleaseContext(m_context);

#ifdef GRALE_DEBUG_OPENCLSENTINEL
	setSavedQueue(nullptr);
#endif // GRALE_DEBUG_OPENCLSENTINEL

	m_context = 0;
	m_queue = 0;
	for (auto &p : m_programs)
		p = 0;
	for (auto &k : m_kernels)
		k = 0;
	m_device = 0;
}

template<int numKernels> 
int OpenCLMultiKernel<numKernels>::getRotatedDeviceIndex()
{
	std::string fileName = "/dev/shm/grale_mpopencl_nextdevice.dat";
	grale::getenv("GRALE_OPENCL_AUTODEVICEFILE", fileName); // Doesn't change file name if envvar not set

	m_perNodeCounter = std::make_unique<grale::PerNodeCounter>(fileName);

	int idx = m_perNodeCounter->getCount();
	if (idx < 0)
	{
		m_perNodeCounter.reset();
		setErrorString("Couldn't read per-node device index from file '" + fileName + "': " + m_perNodeCounter->getErrorString());
		return -1;
	}

	int numDevices = getDeviceCount();
	if (numDevices < 0)
	{
		setErrorString("Error getting device count: " + getErrorString());
		return -1;
	}
	if (numDevices == 0)
	{
		setErrorString("Unexpectedly got zero GPU devices");
		return -1;
	}
	return idx%numDevices;
}
