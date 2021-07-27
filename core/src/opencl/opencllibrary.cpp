#ifdef GRALE_LOADLIBRARY
#include <windows.h>
#else
#include <dlfcn.h>
#endif // GRALE_LOADLIBRARY
#include "opencllibrary.h"
#include "utils.h"

using namespace std;

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
	GETFUNCTION(clGetPlatformInfo)
	GETFUNCTION(clGetDeviceInfo)

	return true;
}

