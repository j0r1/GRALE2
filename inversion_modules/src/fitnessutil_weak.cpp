#include "fitnessutil.h"
#include "fitnesscomponent.h"
#include <grale/projectedimagesinterface.h>
#include <assert.h>
#include <limits>
#include <complex>
#include <iomanip>

using namespace std;

namespace grale
{

inline float calculateEllipticityProbability(float epsilon, float f1, float f2, 
											float axx, float ayy, float axy, float distFrac,
											const DiscreteFunction<float> &baDistFunction)
{
	float gamma1 = 0.5f*(axx-ayy)*distFrac;
	float gamma2 = axy*distFrac;
	float kappa = 0.5f*(axx+ayy)*distFrac;
	
	float oneMinusKappa = 1.0f-kappa;
	if (ABS(oneMinusKappa) < epsilon) // avoid division by zero
			oneMinusKappa = copysign(epsilon, oneMinusKappa);

	float g1 = gamma1/oneMinusKappa;
	float g2 = gamma2/oneMinusKappa;

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
	float baProb = baDistFunction(F);
	float F1 = F+1.0f;
	float prob = 1.0f/(float(CONST_PI) * (F + epsilon) * F1*F1) * baProb * jac; // avoid div by zero
	return prob;
}

float calculateEllipticityProbabilityWithError(float epsilon, float startFromSigmaFactor, float sigmaSteps,
											float f1, float f2, float sigma1, float sigma2,
											float axx, float ayy, float axy, float distFrac,
											const DiscreteFunction<float> &baDistFunction)
{
	if (sigma1 == 0 && sigma2 == 0)
		return calculateEllipticityProbability(epsilon, f1, f2, axx, ayy, axy, distFrac, baDistFunction);

	float f1Start = f1 - sigma1*startFromSigmaFactor;
	float f1End = f1 + sigma1*startFromSigmaFactor;
	float f2Start = f2 - sigma2*startFromSigmaFactor;
	float f2End = f2 + sigma2*startFromSigmaFactor;
	float df1 = (f1End-f1Start)/((float)(sigmaSteps-1));
	float df2 = (f2End-f2Start)/((float)(sigmaSteps-1));
	float weightSum = 0;
	float probSum = 0;

	for (int y = 0 ; y < sigmaSteps ; y++)
	{
		float f2Sample = f2Start + y*df2;
		float difff2 = (f2Sample - f2)/sigma2;
		float gauss2 = std::exp(-0.5f*difff2*difff2);

		for (int x = 0 ; x < sigmaSteps ; x++)
		{
			float f1Sample = f1Start + x*df1;
			float difff1 = (f1Start - f1)/sigma1;
			float gauss1 = std::exp(-0.5f*difff1*difff1);

			float factor = gauss1*gauss2;
			weightSum += factor;

			float prob = 0.0f;
			if (f1Sample*f1Sample + f2Sample*f2Sample < 1.0f)
				prob = calculateEllipticityProbability(epsilon, f1Sample, f2Sample, axx, ayy, axy, distFrac, baDistFunction);

			probSum += factor * prob;
		}
	}
	return probSum/weightSum;
};

// This is based on a Bayesian calculation, assuming the ellipticities
// are known without error. Further assuming a uniform 
// prior on the basis function weights.
// Then, only the jacobians of the ellipSrc to ellipImg transforms
// weighted by possibly unknown distance fractions

// TODO: this fitness can (and will) get negative. The algorithm does
//       seem to keep working, but for now I'm not really sure if there
//       are unintended consequences of this.

// For the bayesian version, we assume that Dds/Ds of the extended images
// was set to one; the actual Dds/Ds should be stored in the shear weights
// entries. If set to zero, the actual redshift is unknown and a weighted
// average will be used based on distanceFractionWeights (assumed to be
// normalized)

void dumpFunction(const DiscreteFunction<float> &f, const string &fname)
{
	float x0, x1;
	f.getLimits(x0, x1);

	ofstream s(fname, ofstream::out);
	if (!s.is_open())
		cerr << "ERROR: couldn't open " << fname << endl;
	
	const int N = 1024;
	for (int i = 0 ; i < N ; i++)
	{
		float factor = (float)i/(float)(N-1);
		float x = (1.0f-factor)*x0 + factor*x1;
		float y = f(x);

		s << std::setprecision(8) << x << " " << y << endl;
	}
}

float calculateWeakLensingFitness_Bayes(const ProjectedImagesInterface &interface,
	const vector<int> &weakIndices,
	const std::vector<int> &avgDensPriorIndices,
	const vector<vector<float>> &preCalcDistFrac,
	const DiscreteFunction<float> &distFracFunction,
	const std::vector<std::pair<float,float>> &zDistDistFracAndProb,
	const DiscreteFunction<float> &baDistFunction,
	float startFromSigmaFactor, int sigmaSteps,
	float zLens)
{
	// // TODO: for debugging
	// static bool first = true;
	// if (first)
	// {
	// 	cerr << "DEBUG: z_lens = " << zLens << endl;
	// 	first = false;
	// 	dumpFunction(distFracFunction, "distfracfunction.txt");
	// 	dumpFunction(baDistFunction, "badistfunction.txt");

	// 	ofstream f("distfractandprob.txt", ofstream::out);
	// 	if (f.is_open())
	// 	{
	// 		for (auto dp : zDistDistFracAndProb)
	// 			f << setprecision(8) << dp.first << " " << dp.second << endl;
	// 	}
	// }

	const float epsilon = 1e-6; // to avoid division by zero
	float shearFitness = 0;

	assert(weakIndices.size() == preCalcDistFrac.size());
	for (int sIdx = 0 ; sIdx < weakIndices.size() ; sIdx++)
	{
		const int s = weakIndices[sIdx];

		assert(s >= 0 && s < interface.getNumberOfSources());
		int numPoints = interface.getNumberOfImagePoints(s);
		const float *pEll1 = interface.getOriginalProperties(ImagesData::ShearComponent1, s);
		const float *pEll2 = interface.getOriginalProperties(ImagesData::ShearComponent2, s);
		const float *pZ = interface.getOriginalProperties(ImagesData::Redshift, s);
		const float *pZSigma = interface.getOriginalProperties(ImagesData::RedshiftUncertainty, s);
		const float *pAxx = interface.getDerivativesXX(s);
		const float *pAyy = interface.getDerivativesYY(s);
		const float *pAxy = interface.getDerivativesXY(s);
		const float *pEllError1 = (interface.hasOriginalProperty(ImagesData::ShearUncertaintyComponent1, s))?interface.getOriginalProperties(ImagesData::ShearUncertaintyComponent1, s):nullptr;
		const float *pEllError2 = (interface.hasOriginalProperty(ImagesData::ShearUncertaintyComponent2, s))?interface.getOriginalProperties(ImagesData::ShearUncertaintyComponent2, s):nullptr;
		
		const float *pDistFrac = preCalcDistFrac[sIdx].data(); // Note that we need sIdx here, not s
		assert(preCalcDistFrac[sIdx].size() == numPoints);

		// TODO: for debugging
		// {
		// 	serut::FileSerializer fser;
		// 	fser.open("baddump.dat", serut::FileSerializer::ReadOnly);
		// 	int32_t pts;
	
		// 	fser.readInt32(&pts);
		// 	fser.readFloats((float*)pEll1, numPoints);
		// 	fser.readFloats((float*)pEll2, numPoints);
		// 	fser.readFloats((float*)pDistFrac, numPoints);
		// 	fser.readFloats((float*)pAxx, numPoints);
		// 	fser.readFloats((float*)pAyy, numPoints);
		// 	fser.readFloats((float*)pAxy, numPoints);
		// 	cerr << "BAD FITNESS SETTINGS LOADED" << endl;
		// }

		assert(pEll1 && pEll2 && pAxx && pAyy && pAxy);
		
		for (int i = 0 ; i < numPoints ; i++)
		{
			float axx = pAxx[i];
			float ayy = pAyy[i];
			float axy = pAxy[i];
			float distFrac = pDistFrac[i]; // Check later if this is initialized
			float z = pZ[i];
			float zSigma = pZSigma[i];

			float f1 = pEll1[i]; // entries represent measured galaxy ellipticities
			float f2 = pEll2[i];
			float sigma1 = (pEllError1 != nullptr)?pEllError1[i]:0.0f;
			float sigma2 = (pEllError2 != nullptr)?pEllError2[i]:0.0f;

			float elliptProb = 0;
			float weightSum = 0;

			if (z != 0 && zSigma == 0) // Redshift is known accurately
			{
				assert(distFrac >= 0 && distFrac < 1);
				assert(distFracFunction(z) == distFrac);

				elliptProb = calculateEllipticityProbabilityWithError(epsilon, startFromSigmaFactor, sigmaSteps,
					f1, f2, sigma1, sigma2, axx, ayy, axy, distFrac, baDistFunction);
				weightSum = 1.0f;
			}
			else if (z != 0 && zSigma != 0) // Redshift is known with uncertainty
			{

				// // Unknown redshift/distance fraction, use weighted average according to some distribution
				// for (auto fracAndProb : unknownDistFracWeightsNormed)
				// 	elliptProb += calculateEllipticityProbabilityWithError(epsilon, startFromSigmaFactor, sigmaSteps,
				// 		f1, f2, sigma1, sigma2, axx, ayy, axy, fracAndProb.first) * fracAndProb.second;

				float zStart = z - zSigma*startFromSigmaFactor;
				float zEnd = z + zSigma*startFromSigmaFactor;
				float dz = (zEnd-zStart)/((float)(sigmaSteps-1));

				for (int i = 0 ; i < sigmaSteps ; i++)
				{
					float zSample = zStart + i*dz;
					float diffz = (zSample - z)/zSigma;
					float zProb = std::exp(-0.5f*diffz*diffz);
					float distFrac = distFracFunction(zSample);

					if (zSample > zLens) // Otherwise probability is zero
					{
						weightSum += zProb;
						elliptProb += calculateEllipticityProbabilityWithError(epsilon, startFromSigmaFactor, sigmaSteps,
							f1, f2, sigma1, sigma2, axx, ayy, axy, distFrac, baDistFunction) * zProb;
					}
				}
			}
			else if (z == 0 && zSigma == 0)
			{
				// No knowledge about redshift whatsoever
				assert(zDistDistFracAndProb.size() > 0);
				for (auto dp : zDistDistFracAndProb)
				{
					float distFrac = dp.first;
					float zProb = dp.second;

					weightSum += zProb;
					elliptProb += calculateEllipticityProbabilityWithError(epsilon, startFromSigmaFactor, sigmaSteps,
				 		f1, f2, sigma1, sigma2, axx, ayy, axy, distFrac, baDistFunction) * zProb;
				}
			}
			else
				assert(0);

			elliptProb /= weightSum;
			if (elliptProb <= 0) // avoid problems with log
				elliptProb = epsilon; 
			
			float logElliptProb = std::log(elliptProb);
			
			// // TODO: for debugging
			// if (isinf(logElliptProb) || isnan(logElliptProb))
			// {
			// 	serut::FileSerializer fser;
			// 	fser.open("baddump.dat", serut::FileSerializer::WriteOnly);
			// 	fser.writeInt32(numPoints);
			// 	fser.writeFloats(pEll1, numPoints);
			// 	fser.writeFloats(pEll2, numPoints);
			// 	fser.writeFloats(pDistFrac, numPoints);
			// 	fser.writeFloats(pAxx, numPoints);
			// 	fser.writeFloats(pAyy, numPoints);
			// 	fser.writeFloats(pAxy, numPoints);

			// 	serut::FileSerializer fser2;
			// 	vector<float> settings { epsilon, f1, f2, axx, ayy, axy, distFrac };
			// 	fser2.open("badsettings.dat", serut::FileSerializer::WriteOnly);
			// 	fser2.writeFloats(settings);
			// 	cerr << "BAD FITNESS DETECTED, BAILING. Point is " << i << endl;
			// 	exit(-1);
			// }
			
			shearFitness += -logElliptProb; // use negative to search for a minimum
		}
	}

	float avgKappa = 0;
	int numKappaPoints = 0;
	for (int s : avgDensPriorIndices)
	{
		int numPoints = interface.getNumberOfImagePoints(s);
		const float *pKappa = interface.getConvergence(s);

		for (int i = 0 ; i < numPoints ; i++)
			avgKappa += pKappa[i];
		numKappaPoints += numPoints;
	}

	avgKappa /= (float)numKappaPoints;
	const int kappaMin = epsilon; // TODO: make this configurable? Avoid log of zero
	if (avgKappa < kappaMin)
		avgKappa = kappaMin;

	// prob(dens) = cte * 1/dens => logprob = - log(dens) + const
	// since we're using the negative log prob, it's just log(dens)
	shearFitness += log(avgKappa);

	// static bool first = true;
	// if (first)
	// {
	// 	first = false;
	// 	cerr << "weakIndices.size = " << weakIndices.size() << endl;
	// 	cerr << "avgDensPriorIndices.size = " << avgDensPriorIndices.size() << endl;
	// 	cerr << "numKappaPoints = " << numKappaPoints << endl;
	// 	cerr << "avgKappa = " << avgKappa << endl;
	// }

	assert(!isnan(shearFitness));
	return shearFitness;
}

float calculateWeakLensingFitness(const ProjectedImagesInterface &interface, const vector<int> &weakIndices,
								  WeakLensingType type, const vector<float> &oneMinusKappaThreshold)
{
	assert(weakIndices.size() == oneMinusKappaThreshold.size());

	const float epsilon = 1e-6; // to avoid division by zero
	float shearFitness = 0;
	int usedPoints = 0;

	for (int sIdx = 0 ; sIdx < weakIndices.size() ; sIdx++)
	{
		const int s = weakIndices[sIdx];

		assert(s >= 0 && s < interface.getNumberOfSources());
		int numPoints = interface.getNumberOfImagePoints(s);
		const float *pStoredShear1 = interface.getOriginalProperties(ImagesData::ShearComponent1, s);
		const float *pStoredShear2 = interface.getOriginalProperties(ImagesData::ShearComponent2, s);
		const float *pWeights = interface.getOriginalProperties(ImagesData::ShearWeight, s);
		const float *pCalcShear1 = interface.getShearComponents1(s);
		const float *pCalcShear2 = interface.getShearComponents2(s);
		const float *pConvergence = interface.getConvergence(s);

		assert(pStoredShear1 && pStoredShear2 && pCalcShear1 && pCalcShear2 && pConvergence && pWeights);
		float threshold = oneMinusKappaThreshold[sIdx];
		
		for (int i = 0 ; i < numPoints ; i++)
		{
			float gamma1 = pCalcShear1[i]; 
			float gamma2 = pCalcShear2[i];
			float kappa = pConvergence[i];
			float oneMinusKappa = 1.0f-kappa;
			float weight = pWeights[i];

			// Note: we're also using the threshold here in case regular shear is used.
			//       This probably doesn't make much sense but I leave it to the user to
			//       specify that the oneMinusKappaThreshold should be zero in that case.
			if (ABS(oneMinusKappa) >= threshold)
			{
				if (ABS(oneMinusKappa) < epsilon) // avoid division by zero
					oneMinusKappa = copysign(epsilon, oneMinusKappa);

				float d1 = 0, d2 = 0;
				if (type == RealShear) // Real shear, for testing
				{
					// In this case, despite the naming, pStoredShear is supposed to
					// hold the 'normal' shear and not the reduced shear
					d1 = (gamma1-pStoredShear1[i]);
					d2 = (gamma2-pStoredShear2[i]);
				}
				else // Reduced
				{
					float g1 = gamma1/oneMinusKappa;
					float g2 = gamma2/oneMinusKappa;

					if (type == RealReducedShear) // use just reduced shear (for testing)
					{
						d1 = (g1-pStoredShear1[i]);
						d2 = (g2-pStoredShear2[i]);

						shearFitness += (d1*d1 + d2*d2)*weight;
					}
					else if (type == AveragedEllipticities) // Ellipticities
					{
						assert(type == AveragedEllipticities);
						complex<float> Ee { pStoredShear1[i], pStoredShear2[i] }; // entry represents averaged ellipticities
						complex<float> g { g1, g2 };

						// See eg https://arxiv.org/pdf/astro-ph/0509252.pdf
						// E(e) = g if    |g| <= 1
						// E(e) = 1/g* if |g| > 1
						
						if (std::abs(g) > 1.0f)
							g = 1.0f/std::conj(g);

						d1 = (g.real()-Ee.real());
						d2 = (g.imag()-Ee.imag());
						
						shearFitness += (d1*d1 + d2*d2)*weight;
					}
					else // BayesianEllipticities is covered somewhere else, signal error
					{
						return numeric_limits<float>::quiet_NaN();
					}
				}

				usedPoints++;
			}
		}
	}

	if (usedPoints == 0)
		shearFitness = 1e30; // Avoid creating a solution that dominates this fitness because there are no points
	else
		shearFitness /= usedPoints; // this also covers all the weak lensing data sets

	assert(!isnan(shearFitness));
	return shearFitness;
}

} // end namespace
