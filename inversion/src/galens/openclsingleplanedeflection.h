#pragma once

#include "graleconfig.h"
#include "vector2d.h"
#include "openclmultikernel.h"
#include "oclutils.h"
#include <eatk/vectorgenomefitness.h>
#include <errut/booltype.h>
#include <vector>
#include <string>
#include <memory>
#include <mutex>
#include <set>
#include <unordered_map>
#include <limits>

namespace grale
{

// For parametric lens inversion

// Idea is to prepare kernel to calculate alphas, (derivatives, potentials later, if needed)
// for a number of parameter sets (number of genomes for example)
// Here we just supply a number of deflection points, filtering unique ones (because of same
// null space grids for example) will be done before
// These deflection points will be scaled float values
class OpenCLSinglePlaneDeflection
{
public:
	OpenCLSinglePlaneDeflection();
	~OpenCLSinglePlaneDeflection();

	errut::bool_t init(const std::vector<Vector2Df> &thetas, // already transformed into the correct units
					   const std::vector<float> &thetaUncert, // may be empty, otherwise in correct units and same length as thetas 
					   const std::vector<int> &templateIntParameters, // these cannot change
					   const std::vector<float> &templateFloatParameters, // only floating point params can change
					   const std::vector<size_t> changeableParameterIndices,
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
					   int devIdx, // negative means rotate
					   uint64_t initialUncertSeed
					   );

	void destroy();
	int getDeviceIndex() const { return m_devIdx; }

	errut::bool_t calculateDeflection(const std::vector<float> &parameters,
									  std::vector<Vector2Df> &allAlphas,
									  std::vector<float> &allAxx,
									  std::vector<float> &allAyy,
									  std::vector<float> &allAxy,
									  std::vector<float> &allPotentials
									  );
	// should have numChangebleParams * numParamSets length
	//     Here we can either modify the full parameters on the CPU and upload
	//     these, or upload only these parameters and let a kernel change them
	//     in the full parameters
protected:
	errut::bool_t randomizeInputPositions();

	bool m_init = false;
	int m_devIdx = -1;
	std::unique_ptr<OpenCLMultiKernel<3>> m_cl;

	oclutils::CLMem m_clThetas;
	oclutils::CLMem m_clIntParams;
	oclutils::CLMem m_clFloatParams; // copy of parameters for each genome
	oclutils::CLMem m_clAllResults;
	oclutils::CLMem m_clChangedParamsBuffer;
	oclutils::CLMem m_clChangeableParamIndices;

	oclutils::CLMem m_clThetaUncerts; // Will be fixed, input parameters
	oclutils::CLMem m_clThetaWithAdditions; // on calculateDeflection, bases on uncert differences for thetas are calculated and stored
	oclutils::CLMem m_clRngStates; // The RNG states that will be used for this

	size_t m_numPoints, m_numFloatParams, m_currentNumParamSets, m_maxNumParamSets;
	std::vector<cl_float> m_floatParamsCopy; // Single float params
	std::vector<cl_float> m_allFloatParams; // repeats of float params
	std::vector<cl_float> m_allResultsBuffer;
	std::vector<size_t> m_changeableParameterIndices;
};

// Using same single instance code as in OpenCLMultiPlaneCalculator (for multiplane)
// TODO: make this common code somehow
class OpenCLSinglePlaneDeflectionInstance : private OpenCLSinglePlaneDeflection
{
public:
	OpenCLSinglePlaneDeflectionInstance();
	~OpenCLSinglePlaneDeflectionInstance();

	static errut::bool_t initInstance(uint64_t userId,const std::vector<Vector2Df> &thetas, // already transformed into the correct units
					   const std::vector<float> &thetaUncert, // may be empty, otherwise in correct units and same length as thetas 
					   const std::vector<int> &templateIntParameters, // these cannot change
					   const std::vector<float> &templateFloatParameters, // only floating point params can change
					   const std::vector<size_t> changeableParameterIndices,
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
					   int devIdx,
					   uint64_t initialUncertSeed
					   );
	static void releaseInstance(uint64_t userId);
	static OpenCLSinglePlaneDeflectionInstance &instance();

	void setTotalGenomesToCalculate(size_t iteration, size_t num);
	errut::bool_t scheduleCalculation(const eatk::FloatVectorGenome &genome);
	bool getResultsForGenome(const eatk::FloatVectorGenome &genome,
	                         std::vector<Vector2Df> &alphas, std::vector<float> &axx,
							 std::vector<float> &ayy, std::vector<float> &axy,
							 std::vector<float> &potential);

	int getRequestedDeviceIndex() const { return m_requestedDevIdx; }
private:
	static std::unique_ptr<OpenCLSinglePlaneDeflectionInstance> s_instance;
	static std::mutex s_instanceMutex;
	static std::set<uint64_t> s_users;
	static bool s_initTried;

	bool m_calculationDone = false;
	size_t m_totalGenomesToCalculate = 0;
	size_t m_prevIteration = std::numeric_limits<size_t>::max();
	// Map a genome pointer to an offset
	std::unordered_map<const eatk::FloatVectorGenome *, size_t> m_genomeOffsets;
	std::mutex m_mutex;
	std::vector<float> m_floatBuffer;

	std::vector<Vector2Df> m_allAlphas;
	std::vector<float> m_allAxx, m_allAyy, m_allAxy;
	std::vector<float> m_allPotentials;
	int m_requestedDevIdx = -1;
};	
}
