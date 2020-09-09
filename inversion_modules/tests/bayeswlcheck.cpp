#include <serut/fileserializer.h>
#include <grale/vector2d.h>
#include <grale/constants.h>
#include <vector>
#include <iostream>
#include <complex>
#include <cmath>

using namespace std;
using namespace grale;

int main(int argc, char const *argv[])
{
	const float epsilon = 1e-6; // to avoid division by zero

	auto calculateEllipticityProbability = [epsilon](float f1, float f2, 
	                                          float axx, float ayy, float axy, float distFrac) -> float
	{
		float gamma1 = 0.5f*(axx-ayy)*distFrac;
		float gamma2 = axy*distFrac;
		float kappa = 0.5f*(axx+ayy)*distFrac;
		
		float oneMinusKappa = 1.0f-kappa;
		if (ABS(oneMinusKappa) < epsilon) // avoid division by zero
				oneMinusKappa = copysign(epsilon, oneMinusKappa);

		float g1 = gamma1/oneMinusKappa;
		float g2 = gamma2/oneMinusKappa;

		// This is based on a Bayesian calculation, assuming the ellipticities
		// are known without error. Further assuming a uniform distribution 
		// for the b/a ratio of the elliptic source shapes, and a uniform 
		// prior on the basis function weights.
		// Then, only the jacobians of the ellipSrc to ellipImg transforms
		// weighted by possibly unknown distance fractions

		// TODO: this fitness can (and will) get negative. The algorithm does
		//       seem to keep working, but for now I'm not really sure if there
		//       are unintended consequences of this.
		
		// weights are ignored currently
		float gSq = g1*g1 + g2*g2;
		float fSq = f1*f1 + f2*f2;
		complex<float> eImg = { f1, f2 };
		complex<float> g = { g1, g2 };
		float F = 0;
		float jacRoot = 0;
		
		// A factor 2 has been omitted, and logprob sign reversed so that we'll be
		// looking for a minimum
		if (gSq <= 1)
		{
			F = abs(eImg - g)/(abs(1.0f - conj(g)*eImg) + epsilon);
			jacRoot = (gSq-1.0f)/(abs(fSq*gSq - 2.0f*f1*g1 - 2.0f*f2*g2 + 1.0f) + epsilon);
		}
		else
		{
			F = abs(1.0f - g*conj(eImg))/(abs(conj(eImg) - conj(g)) + epsilon);
			jacRoot = (gSq-1.0f)/(abs(fSq + gSq - 2.0f*f1*g2 - 2.0f*f2*g2) + epsilon);
		}

		float jac = jacRoot*jacRoot;
        cerr << "jac = " << jac << endl;
		// Note that this is actualy still a proportionality, depending on the b/a distribution
		// With this proportionality, every b/a, from 0 to 1 is possible with equal probability
		float F1 = F+1.0f;
        cerr << "F = " << F << endl;
        cerr << "F1 = " << F1 << endl;
		float prob = 1.0f/(float(CONST_PI) * (F + epsilon) * F1*F1)*jac; // avoid div by zero
		return prob;
	};

    serut::FileSerializer fser2;
    vector<float> settings(7);
    if (!fser2.open("badsettings.dat", serut::FileSerializer::ReadOnly))
        cerr << "Coudln't open file: " << fser2.getErrorString() << endl;

    if (!fser2.readFloats(settings))
        cerr << "Coudln't read floats: " << fser2.getErrorString() << endl;

    float fitness = calculateEllipticityProbability(settings[1], settings[2], settings[3], settings[4],
                                                    settings[5], settings[6]);

    cerr << fitness << " " << (-std::log(fitness)) << endl;
    return 0;
}