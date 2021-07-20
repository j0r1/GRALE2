#include "openclkernel.h"
#include <iostream>
#include <vector>

using namespace std;

class OCLTest : public OpenCLKernel
{
public:
	OCLTest() { }
	~OCLTest() { }

	bool run();
};

bool OCLTest::run()
{
	cl_uint numPlatforms;
	cl_int err = clGetPlatformIDs(0, nullptr, &numPlatforms);
	if (err != CL_SUCCESS)
	{
		setErrorString("Can't get platform IDs:" + getCLErrorString(err));
		return false;
	}
	cout << "Detected " << numPlatforms << " platforms" << endl;
	if (numPlatforms == 0)
	{
		setErrorString("No platforms available");
		return false;
	}

	vector<cl_platform_id> platforms(numPlatforms);
	clGetPlatformIDs(numPlatforms, platforms.data(), nullptr);
	cout << "Platforms: " << endl;
	for (auto x : platforms)
	{
		char name[1024] = { 0 };
		err = clGetPlatformInfo(x, CL_PLATFORM_NAME, sizeof(name)-1, name, nullptr);
		if (err != CL_SUCCESS)
		{
			setErrorString("Can't get platform info: " + getCLErrorString(err));
			return false;
		}
		cout << " " << x << ": " << name << endl;
	}

	cl_platform_id platform = platforms[0];
	cout << "Using first platform " << platform << endl;

	cl_uint numDevices = 0;
	err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &numDevices);
	if (err != CL_SUCCESS || numDevices == 0)
	{
		setErrorString("Can't find any GPU devices");
		return false;
	}

	cout << "Found " << numDevices << " GPU devices" << endl;
	vector<cl_device_id> devices(numDevices);
	clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, numDevices, devices.data(), nullptr);
	cout << "GPU Devices:" << endl;
	for (auto x : devices)
	{
		char name[1024] = { 0 };
		err = clGetDeviceInfo(x, CL_DEVICE_NAME, sizeof(name)-1, name, nullptr);
		if (err != CL_SUCCESS)
		{
			setErrorString("Can't get device info: " + getCLErrorString(err));
			return false;
		}
		cout << " " << x << ": " << name << endl;
	}

	return true;
}

int main(int argc, char *argv[])
{
	string lib = "/usr/lib64/libOpenCL.so";
	if (argc > 1)
		lib = string(argv[1]);

	cout << "Using library " << lib << endl;

	OCLTest tst;
	if (!tst.loadLibrary(lib))
	{
		cerr << "Can't load library: " << tst.getErrorString() << endl;
		return -1;
	}

	if (!tst.run())
	{
		cerr << "Error running test: " << tst.getErrorString() << endl;
		return -1;
	}
	return 0;
}
