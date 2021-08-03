#ifdef GRALE_LOADLIBRARY
#include <windows.h>
#else
#include <dlfcn.h>
#endif // GRALE_LOADLIBRARY
#include "opencllibrary.h"
#include "utils.h"
#include <errut/booltype.h>
#include <sstream>
#include <vector>

using namespace std;
using namespace errut;

#ifdef GRALE_LOADLIBRARY
#define LOADLIBRARY(x) (void*)LoadLibrary(x)
#define SYMBOL(l, x) GetProcAddress((HMODULE)l, x)
#define CLOSELIBRARY(x) FreeLibrary((HMODULE)x)
#define GETERRORSTRING() ("Error code " + to_string((int)GetLastError()))
#else
#define LOADLIBRARY(x) dlopen(x, RTLD_LAZY)
#define SYMBOL(l, x) dlsym(l, x)
#define CLOSELIBRARY(x) dlclose(x)
#define GETERRORSTRING() (string(dlerror()))
#endif // GRALE_LOADLIBRARY

string OpenCLLibrary::getLibraryName()
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
    grale::getenv("GRALE_OPENCLLIB", library); // Doesn't change the value if environment variable not set

	return library;
}

OpenCLLibrary::OpenCLLibrary()
{
	clBuildProgram = nullptr;
	clCreateCommandQueue = nullptr;
	clCreateContext = nullptr;
	clCreateKernel = nullptr;
	clCreateProgramWithSource = nullptr;
	clGetDeviceIDs = nullptr;
	clGetPlatformIDs = nullptr;
	clGetPlatformInfo = nullptr;
	clGetDeviceInfo = nullptr;
	clGetProgramBuildInfo = nullptr;
	clReleaseCommandQueue = nullptr;
	clReleaseContext = nullptr;
	clReleaseKernel = nullptr;
	clReleaseProgram = nullptr;
	clCreateBuffer = nullptr;
	clSetKernelArg = nullptr;
	clReleaseMemObject = nullptr;
	clEnqueueNDRangeKernel = nullptr;
	clFinish = nullptr;
	clEnqueueReadBuffer = nullptr;
	clEnqueueWriteBuffer = nullptr;
	clSetEventCallback = nullptr;
	clGetEventInfo = nullptr;
	clFlush = nullptr;

	m_pModule = nullptr;
}

OpenCLLibrary::~OpenCLLibrary()
{
	if (m_pModule)
	{
		CLOSELIBRARY(m_pModule);
		m_pModule = nullptr;
	}
}

bool OpenCLLibrary::loadLibrary(const std::string &libraryName)
{
	if (m_pModule)
	{
		setErrorString("A library has already been opened");
		return false;
	}

	m_pModule = LOADLIBRARY(libraryName.c_str());
	if (m_pModule == nullptr)
	{
		setErrorString("Unable to open OpenCL library: " + GETERRORSTRING());
		return false;
	}

#define GETFUNCTION(x) \
	x = (decltype(x))(SYMBOL(m_pModule, #x)); \
	if (x == nullptr) \
	{ \
		setErrorString("Unable to locate symbol " #x " in OpenCL library"); \
		CLOSELIBRARY(m_pModule); \
		m_pModule = nullptr; \
		return false; \
	}

	GETFUNCTION(clBuildProgram)
	GETFUNCTION(clCreateCommandQueue)
	GETFUNCTION(clCreateContext)
	GETFUNCTION(clCreateKernel)
	GETFUNCTION(clCreateProgramWithSource)
	GETFUNCTION(clGetDeviceIDs)
	GETFUNCTION(clGetPlatformIDs)
	GETFUNCTION(clGetProgramBuildInfo)
	GETFUNCTION(clReleaseCommandQueue)
	GETFUNCTION(clReleaseContext)
	GETFUNCTION(clReleaseKernel)
	GETFUNCTION(clReleaseProgram)
	GETFUNCTION(clCreateBuffer)
	GETFUNCTION(clSetKernelArg)
	GETFUNCTION(clReleaseMemObject)
	GETFUNCTION(clEnqueueNDRangeKernel)
	GETFUNCTION(clFinish)
	GETFUNCTION(clEnqueueReadBuffer)
	GETFUNCTION(clEnqueueWriteBuffer)
	GETFUNCTION(clGetPlatformInfo)
	GETFUNCTION(clGetDeviceInfo)
	GETFUNCTION(clSetEventCallback)
	GETFUNCTION(clGetEventInfo)
	GETFUNCTION(clFlush)

	return true;
}

bool OpenCLLibrary::getPlatformAndDeviceCount(cl_platform_id &platformId, int &deviceCount) const
{
	if (!isOpen())
	{
		setErrorString("No OpenCL library has been initialized yet");
		return false;
	}

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

	platformId = platform;
	deviceCount = (int)numDevices;
	return true;
}

int OpenCLLibrary::getDeviceCount() const
{
	cl_platform_id platform;
	int count;
	if (!getPlatformAndDeviceCount(platform, count))
		return -1;
	return count;
}

string OpenCLLibrary::getCLErrorString(int errNum)
{
	return "OpenCL error code " + to_string(errNum);
}

