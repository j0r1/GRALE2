#pragma once

#include "graleconfig.h"

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

	// bool_t init(const vector<Vector2Df> &thetas,
	//             const vector<float> &templateParameters,
	//             const vector<size_t> changeableParameterIndices,
	//             size_t numParamSets, // number of genomes for example
	//             const string &deflectionCode
	//             )
	//
	// bool_t calculateDeflection(const vector<float> &parameters) // should have numChangebleParams * numParamSets length
	//     Here we can either modify the full parameters on the CPU and upload
	//     these, or upload only these parameters and let a kernel change them
	//     in the full parameters
};

}
