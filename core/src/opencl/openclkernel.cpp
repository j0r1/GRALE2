#include "openclkernel.h"
#include <stdio.h>
#include <iostream>

using namespace std;

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

bool OpenCLKernel::init()
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

	cl_platform_id platforms[16];
	cl_uint numPlatforms, numDevices;
	cl_int err = clGetPlatformIDs(16, platforms, &numPlatforms);

	if (err != CL_SUCCESS)
	{
		setErrorString(getCLErrorString(err));
		releaseAll();
		return false;
	}

	if (numPlatforms < 1)
	{
		setErrorString("No OpenCL platforms available");
		releaseAll();
		return false;
	}

	cl_platform_id platform = platforms[0]; // TODO: make this configurable?

	err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
	if (err != CL_SUCCESS)
	{
		setErrorString("Can't find any GPU devices");
		releaseAll();
		return false;
	}

	m_pDevices = make_unique<cl_device_id[]>(numDevices);

	clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, m_pDevices.get(), NULL);

	m_context = clCreateContext(0, 1, m_pDevices.get(), NULL, NULL, &err);
	if (err != CL_SUCCESS)
	{
		setErrorString(string("Can't create OpenCL context: ") + getCLErrorString(err));
		releaseAll();
		return false;
	}

	m_deviceIndex = 0; // TODO: make this configurable?

	m_queue = clCreateCommandQueue(m_context, m_pDevices[m_deviceIndex], 0, &err);
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

		err = clBuildProgram(m_program, 1, &(m_pDevices[m_deviceIndex]), "", 0, 0);
		if (err != CL_SUCCESS)
		{
			size_t logLength;

			clGetProgramBuildInfo(m_program, m_pDevices[m_deviceIndex], CL_PROGRAM_BUILD_LOG, 0, 0, &logLength);

			unique_ptr<char []> pLog = make_unique<char[]>(logLength+1);

			clGetProgramBuildInfo(m_program, m_pDevices[m_deviceIndex], CL_PROGRAM_BUILD_LOG, logLength, pLog.get(), 0);
			pLog[logLength] = 0;

			failLog = string(pLog.get());

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
	m_pDevices = nullptr;
}

string OpenCLKernel::getCLErrorString(int errNum)
{
	return "OpenCL error code " + to_string(errNum);
}

