#include "openclkernel.h"
#include "utils.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
using namespace errut;

OpenCLKernel::OpenCLKernel()
{
	m_init = false;
	m_context = 0;
	m_queue = 0;
	m_program = 0;
	m_kernel = 0;
}

OpenCLKernel::~OpenCLKernel()
{
	destroy();
}

bool OpenCLKernel::init(int devIdx)
{
	if (!isOpen())
	{
		setErrorString("No OpenCL library has been initialized yet");
		return false;
	}

	if (m_init)
	{
		setErrorString("OpenCL kernel loader is already initialized");
		return false;
	}

	m_currentProgram = "";
	m_currentKernel = "";

	cl_uint numPlatforms;
	cl_int err = clGetPlatformIDs(0, nullptr, &numPlatforms);
	if (err != CL_SUCCESS)
	{
		setErrorString("Can't get number of platforms:" + getCLErrorString(err));
		return false;
	}

	if (numPlatforms == 0)
	{
		setErrorString("No platforms available");
		return false;
	}

	vector<cl_platform_id> platforms(numPlatforms);
	cl_platform_id platform;

	clGetPlatformIDs(numPlatforms, platforms.data(), nullptr);
	
	if (numPlatforms == 1)
		platform = platforms[0];
	else
	{
		int platformIdx;
		bool_t r = grale::getenv("GRALE_OPENCL_PLATFORM", platformIdx, 0, (int)numPlatforms-1);

		if (!r)
		{
			stringstream ss;

			ss << "More than one (" << numPlatforms << ") OpenCL platforms detected, can't read GRALE_OPENCL_PLATFORM environment variable for the platform index: ";
			ss << r.getErrorString() << ". Available platforms are";

			for (size_t i = 0 ; i < platforms.size() ; i++)
			{
				ss << " [" << i << "] ";

				char name[1024] = { 0 };
				err = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(name)-1, name, nullptr);
				if (err != CL_SUCCESS)
					ss << "(unknown)";
				ss << name;
			}

			setErrorString(ss.str());
			return false;
		}
		platform = platforms[platformIdx];
	}

	cl_uint numDevices = 0;
	err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &numDevices);
	if (err != CL_SUCCESS || numDevices == 0)
	{
		setErrorString("Can't find any GPU devices");
		return false;
	}

	if (devIdx < 0 || devIdx >= (int)numDevices)
	{
		setErrorString("Invalid device index (" + to_string(devIdx) + "), there are " + to_string(numDevices) + " detected");
		return false;
	}

	vector<cl_device_id> devices(numDevices);
	clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices.data(), nullptr);
	
	m_device = devices[devIdx];
	m_context = clCreateContext(nullptr, 1, &m_device, nullptr, nullptr, &err);
	if (err != CL_SUCCESS)
	{
		setErrorString(string("Can't create OpenCL context: ") + getCLErrorString(err));
		releaseAll();
		return false;
	}

	m_queue = clCreateCommandQueue(m_context, m_device, 0, &err);
	if (err != CL_SUCCESS)
	{
		setErrorString(string("Can't create OpenCL command queue: ") + getCLErrorString(err));
		releaseAll();
		return false;
	}

	m_init = true;

	return true;
}

bool OpenCLKernel::loadKernel(const string &programString, const string &kernelName, string &failLog)
{
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

	if (m_program == 0 || (m_program != 0 && programString != m_currentProgram))
	{
		// compile
		
		if (m_program)
			clReleaseProgram(m_program);
		if (m_kernel)
			clReleaseKernel(m_kernel);

		m_program = 0;
		m_kernel = 0;
		m_currentProgram = "";
		m_currentKernel = "";

		const char *pProgStr = programString.c_str();
		cl_int err;

		m_program = clCreateProgramWithSource(m_context, 1, &pProgStr, 0, &err);
		if (err != CL_SUCCESS)
		{
			setErrorString(string("Unable to create OpenCL program from specified source: ") + getCLErrorString(err));
			return false;
		}

		err = clBuildProgram(m_program, 1, &m_device, "", nullptr, nullptr);
		if (err != CL_SUCCESS)
		{
			size_t logLength;

			clGetProgramBuildInfo(m_program, m_device, CL_PROGRAM_BUILD_LOG, 0, 0, &logLength);

			vector<char> pLog(logLength+1);
			clGetProgramBuildInfo(m_program, m_device, CL_PROGRAM_BUILD_LOG, logLength, pLog.data(), 0);
			pLog[logLength] = 0;

			failLog = string(pLog.data());

			clReleaseProgram(m_program);
			m_program = 0;

			setErrorString(string("Unable to build OpenCL program: ") + getCLErrorString(err));

			return false;
		}

		m_currentProgram = programString;
	}

	if (m_kernel == 0 || (m_kernel != 0 && m_currentKernel != kernelName))
	{
		cl_int err;

		if (m_kernel)
			clReleaseKernel(m_kernel);

		m_currentKernel = "";

		m_kernel = clCreateKernel(m_program, kernelName.c_str(), &err);
		if (err != CL_SUCCESS)
		{
			setErrorString(string("Unable to get OpenCL kernel from program: ") + getCLErrorString(err));
			return false;
		}

		m_currentKernel = kernelName;
	}

	return true;
}

void OpenCLKernel::destroy()
{
	m_init = false;

	m_currentProgram = "";
	m_currentKernel = "";

	releaseAll();
}

void OpenCLKernel::releaseAll()
{
	if (m_kernel)
		clReleaseKernel(m_kernel);
	if (m_program)
		clReleaseProgram(m_program);
	if (m_queue)
		clReleaseCommandQueue(m_queue);
	if (m_context)
		clReleaseContext(m_context);

	m_context = 0;
	m_queue = 0;
	m_program = 0;
	m_kernel = 0;
	m_device = 0;
}

string OpenCLKernel::getCLErrorString(int errNum)
{
	return "OpenCL error code " + to_string(errNum);
}

