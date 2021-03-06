#ifndef WIN32
#include "multiplanecuda.h"
#include "pernodecounter.h"
#include "utils.h"
#include <dlfcn.h>
#include <sys/time.h>
#include <sstream>
#include <iostream>

using namespace std;

namespace grale
{

MultiPlaneCUDA::MultiPlaneCUDA()
{
	zero();
}

void MultiPlaneCUDA::zero()
{
	mpcuInitMultiPlaneCalculation = nullptr;
	mpcuCalculateSourcePositions = nullptr;
	mpcuGetSourcePositions = nullptr;
	mpcuClearContext = nullptr;
	m_pLibrary = nullptr;
	m_pContext = nullptr;
	m_perNodeCounter.reset();
}

void MultiPlaneCUDA::cleanup()
{
	if (m_pContext)
		mpcuClearContext(m_pContext);
	if (m_pLibrary)
		dlclose(m_pLibrary);
	m_perNodeCounter.reset();
	zero();
}

MultiPlaneCUDA::~MultiPlaneCUDA()
{
	cleanup();
}

#define GETFUNCTION(x) do { \
		x = (decltype(x))dlsym(m_pLibrary, #x ); \
	    if (!x) \
		{ \
			setErrorString("Unable to get library symbol '" #x "':" + string(dlerror())); \
			cleanup(); \
			return false; \
		} \
    } while(0)

bool MultiPlaneCUDA::init(const std::string &libraryPath, int devIdx,
	double angularUnit,
	double h, double W_m, double W_r, double W_v, double w,
	const std::vector<float> &lensRedshifts,
	const std::vector<std::vector<PlummerInfo>> &fixedPlummerParameters, 
	const std::vector<float> &sourceRedshifts,
	const std::vector<std::vector<Vector2Df>> &theta)
{
	// The helper library that is used appears to need the lens redshifts to
	// be provided in an ordered fashion, otherwise incorrect results are
	// generated
	for (size_t i = 1 ; i < lensRedshifts.size() ; i++)
	{
		if (lensRedshifts[i] <= lensRedshifts[i-1])
		{
			setErrorString("The lenses need to be provided with strictly increasing redshifts");
			return false;
		}
	}

	if (m_pLibrary)
	{
		setErrorString("Already initialized");
		return false;
	}

	m_pLibrary = dlopen(libraryPath.c_str(), RTLD_NOW);
	if (!m_pLibrary)
	{
		setErrorString("Unable to open library '" + libraryPath + "': " + string(dlerror()));
		return false;
	}

	GETFUNCTION(mpcuInitMultiPlaneCalculation);
	GETFUNCTION(mpcuCalculateSourcePositions);
	GETFUNCTION(mpcuGetSourcePositions);
	GETFUNCTION(mpcuClearContext);
	GETFUNCTION(mpcuGetNumberOfDevices);

	int numDevices = mpcuGetNumberOfDevices();
	if (numDevices <= 0)
	{
		setErrorString("No CUDA devices appear to be available");
		return false;
	}

	if (devIdx < 0)
	{
		string fileName = "/dev/shm/grale_mpcuda_nextdevice.dat";
		getenv("GRALE_MPCUDA_AUTODEVICEFILE", fileName); // Doesn't change file name if envvar not set

		m_perNodeCounter = make_unique<PerNodeCounter>(fileName);

		int idx = m_perNodeCounter->getCount();
		if (idx < 0)
		{
			setErrorString("Couldn't read per-node device index from file '" + fileName + "': " + m_perNodeCounter->getErrorString());
			m_perNodeCounter.reset();
			return false;
		}

		auto GetTimeStamp = []() {
			struct timeval tv;
			gettimeofday(&tv,NULL);
			return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
		};
		devIdx = idx%numDevices;
		cerr << "DEBUG (" << GetTimeStamp() << "): Automatically using new device index " << devIdx << endl;
	}
	else if (devIdx >= numDevices)
	{
		setErrorString("Invalid device index " + to_string(devIdx) + ", must be in 0..." + to_string(numDevices-1));
		return false;
	}

	int err = mpcuInitMultiPlaneCalculation(angularUnit, h, W_m, W_r, W_v, w,
	                                        lensRedshifts, fixedPlummerParameters, sourceRedshifts, theta,
	                                        &m_pContext, devIdx);
	if (err)
	{
		stringstream ss;
		ss << "Unable to initialize multi plane calculation, got error code: " << err;
		setErrorString(ss.str());
		cleanup();
		return false;
	}

	m_deviceIndex = devIdx; // save the real device index
	m_numSourcePlanes = (int)sourceRedshifts.size();

	return true;
}

bool MultiPlaneCUDA::calculateSourcePositions(const std::vector<std::vector<float>> &massFactors,
                                              const std::vector<float> &sheetDensities)
{
	if (!m_pContext)
	{
		setErrorString("Not initialized");
		return false;
	}
	int err = mpcuCalculateSourcePositions(m_pContext, massFactors, sheetDensities);
	if (err)
	{
		stringstream ss;
		ss << "Got error code: " << err;
		setErrorString(ss.str());
		return false;
	}
	return true;
}

const std::vector<Vector2Df> *MultiPlaneCUDA::getSourcePositions(int srcIdx)
{
	if (!m_pContext)
	{
		setErrorString("Not initialized");
		return nullptr;
	}
	if (srcIdx < 0 || srcIdx >= m_numSourcePlanes)
	{
		setErrorString("Invalid source index");
		return nullptr;
	}
	return &(mpcuGetSourcePositions(m_pContext, srcIdx));
}

} // end namespace

#else 

#include "multiplanecuda.h"
#include "pernodecounter.h"

namespace grale
{

MultiPlaneCUDA::MultiPlaneCUDA()
{
}

void MultiPlaneCUDA::zero()
{
}

void MultiPlaneCUDA::cleanup()
{
}

MultiPlaneCUDA::~MultiPlaneCUDA()
{
}

bool MultiPlaneCUDA::init(const std::string &libraryPath, int devIdx,
	double angularUnit,
	double h, double W_m, double W_r, double W_v, double w,
	const std::vector<float> &lensRedshifts,
	const std::vector<std::vector<PlummerInfo>> &fixedPlummerParameters, 
	const std::vector<float> &sourceRedshifts,
	const std::vector<std::vector<Vector2Df>> &theta)
{
	setErrorString("No Win32 version available");
	return false;
}

bool MultiPlaneCUDA::calculateSourcePositions(const std::vector<std::vector<float>> &massFactors,
	const std::vector<float> &sheetDensities)
{
	setErrorString("No Win32 version available");
	return false;
}

const std::vector<Vector2Df> *MultiPlaneCUDA::getSourcePositions(int srcIdx)
{
	setErrorString("No Win32 version available");
	return nullptr;
}

} // end namespace

#endif // !WIN32
