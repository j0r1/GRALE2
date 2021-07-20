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
	cl_uint numPlatforms, numDevices;
	cl_int err = clGetPlatformIDs(0, nullptr, &numPlatforms);
	if (err != CL_SUCCESS)
	{
		setErrorString("Can't get platform IDs:" + getCLErrorString(err));
		return false;
	}
	cout << "Detected " << numPlatforms << " platforms" << endl;

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
