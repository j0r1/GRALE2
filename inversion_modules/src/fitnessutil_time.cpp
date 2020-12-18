#include "fitnessutil.h"
#include "fitnesscomponent.h"
#include <grale/projectedimagesinterface.h>
#include <assert.h>

using namespace std;

namespace grale
{

float calculateTimeDelayFitnessPaper2009(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices,
		                        const vector<float> &tdScaleFactors)
{
	assert(tdScaleFactors.size() == 0 || tdScaleFactors.size() == sourceIndices.size());

	float timeDelayFitness = 0;

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const float scaleFactor = (tdScaleFactors.size() == 0)?0:tdScaleFactors[sIdx];

		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);
		assert(scaleFactor >= 0);

		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
		float tdfit = 0;
		int count = 0;

		for (int i = 0 ; i < numTimeDelays ; i++)
		{
			int img1, point1;
			float originalDelay1;

			iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

			if (originalDelay1 < -0.0001f) // negative values can be used to simply indicate a source position
				continue;

			for (int j = 0 ; j < numTimeDelays ; j++)
			{
				if (j == i)
					continue;

				int img2, point2;
				float originalDelay2;

				iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

				if (originalDelay2 < -0.0001f)
					continue;

				float realTimeDelayDifference = originalDelay2 - originalDelay1;

				// We have a time delay difference. Now, we'll see how well we can fit this
				// for each backprojected image position

				float squaredSum = 0;

				for (int k = 0 ; k < numTimeDelays ; k++)
				{
					int srcImg, srcPoint;
					float dummyValue;
					Vector2D<float> beta1;

					iface.getOriginalTimeDelay(s, k, &srcImg, &srcPoint, &dummyValue);
					beta1 = iface.getBetas(s, srcImg)[srcPoint];

					for (int l = 0 ; l < numTimeDelays ; l++)
					{
						Vector2D<float> beta2;

						iface.getOriginalTimeDelay(s, l, &srcImg, &srcPoint, &dummyValue);
						beta2 = iface.getBetas(s, srcImg)[srcPoint];

						float delay1 = iface.getTimeDelay(s, img1, point1, beta1);
						float delay2 = iface.getTimeDelay(s, img2, point2, beta2);
						float calculatedDifference = delay2 - delay1;

						float denom = (scaleFactor == 0)?realTimeDelayDifference:scaleFactor;
						float relativeDiff = (realTimeDelayDifference - calculatedDifference)/denom;

						squaredSum += relativeDiff*relativeDiff;
					}
				}

				squaredSum /= (float)(numTimeDelays*numTimeDelays);

				tdfit += squaredSum;
				count++;
			}
		}
		
		if (count)
			tdfit /= (int)count;

		timeDelayFitness += tdfit;
	}

	assert(!isnan(timeDelayFitness));
	return timeDelayFitness;
}

float calculateTimeDelayFitnessNoSrc(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices,
                                             const vector<float> &tdScaleFactors)
{
	assert(tdScaleFactors.size() == 0 || tdScaleFactors.size() == sourceIndices.size());

	float timeDelayFitness = 0;

	double D_d = iface.getLensDistance();
	double z_d = iface.getLensRedshift();
	double dFactor = ((D_d*(1.0+z_d)/SPEED_C) * iface.getAngularScale()*iface.getAngularScale())/(60*60*24);
	float baseFactor = (float)dFactor;

	for (int sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const float scaleFactor = (tdScaleFactors.size() == 0)?0:tdScaleFactors[sIdx];

		assert(iface.getOriginalNumberOfTimeDelays(s) > 0);
		assert(scaleFactor >= 0);

		int numTimeDelays = iface.getOriginalNumberOfTimeDelays(s);
		float tdfit = 0;
		int count = 0;
	
		float factor = baseFactor/iface.getDistanceFraction(s);

		for (int i = 0 ; i < numTimeDelays ; i++)
		{
			int img1, point1;
			float originalDelay1;

			iface.getOriginalTimeDelay(s, i, &img1, &point1, &originalDelay1);

			if (originalDelay1 < -0.0001f) // negative values can be used to simply indicate a source position
				continue;

			Vector2Df alpha1 = iface.getAlphas(s, img1)[point1];
			Vector2Df theta1 = iface.getThetas(s, img1)[point1];
			float phi1 = iface.getLensPotential(s, img1)[point1];

			for (int j = i+1 ; j < numTimeDelays ; j++)
			{
				int img2, point2;
				float originalDelay2;

				iface.getOriginalTimeDelay(s, j, &img2, &point2, &originalDelay2);

				if (originalDelay2 < -0.0001f)
					continue;

				float realTimeDelayDifference = originalDelay2 - originalDelay1;

				Vector2Df alpha2 = iface.getAlphas(s, img2)[point2];
				Vector2Df theta2 = iface.getThetas(s, img2)[point2];
				float phi2 = iface.getLensPotential(s, img2)[point2];

				float phiDiff = (phi2-phi1);
				//float alphaDiff = 0.5*(alpha2.getLengthSquared() - alpha1.getLengthSquared());
				Vector2Df alphaSum = alpha2;
				alphaSum += alpha1;

				Vector2Df thetaDiff = theta2;
				thetaDiff -= theta1;

				float calculatedDifference = (
						0.5*(
							(thetaDiff.getX()*alphaSum.getX()) + (thetaDiff.getY()*alphaSum.getY())
						) 
						-
						phiDiff)*factor;

				float denom = (scaleFactor == 0)?realTimeDelayDifference:scaleFactor;
				float relativeDiff = (realTimeDelayDifference - calculatedDifference)/denom;
					
				tdfit += relativeDiff*relativeDiff;
				count++;
			}
		}
		
		if (count)
			tdfit /= (int)count;

		timeDelayFitness += tdfit;
	}

	assert(!isnan(timeDelayFitness));
	return timeDelayFitness;
}

} // end namespace
