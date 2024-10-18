#pragma once

#include "graleconfig.h"
#include "vector2d.h"
#include "openclmultikernel.h"
#include "oclutils.h"
#include <errut/booltype.h>
#include <vector>
#include <string>
#include <memory>

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
					   const std::vector<int> &templateIntParameters, // these cannot change
					   const std::vector<float> &templateFloatParameters, // only floating point params can change
					   const std::vector<size_t> changeableParameterIndices,
					   //size_t numParamSets, // number of genomes for example
					   const std::string &deflectionKernelCode, const std::string &lensRoutineName,
					   bool uploadFullParameters,
					   int devIdx = 0 // negative means rotate
					   ); // TODO: calculate betas from this as well?

	void destroy();

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
private:
	bool m_init = false;
	bool m_uploadFullParameters;
	std::unique_ptr<OpenCLMultiKernel<3>> m_cl; // TODO how many kernels will we need?

	oclutils::CLMem m_clThetas;
	oclutils::CLMem m_clIntParams;
	oclutils::CLMem m_clFloatParams; // copy of parameters for each genome
	oclutils::CLMem m_clAllResults;
	oclutils::CLMem m_clChangedParamsBuffer;
	oclutils::CLMem m_clChangeableParamIndices;

	size_t m_numPoints, m_numFloatParams, m_currentNumParamSets, m_maxNumParamSets;
	std::vector<cl_float> m_floatParamsCopy; // Single float params
	std::vector<cl_float> m_allFloatParams; // repeats of float params
	std::vector<cl_float> m_allResultsBuffer;
	std::vector<size_t> m_changeableParameterIndices;
};

}
