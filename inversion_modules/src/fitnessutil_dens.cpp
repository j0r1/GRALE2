#include "fitnessutil.h"
#include "pointgroupstorage.h"
#include "fitnesscomponent.h"
#include <grale/projectedimagesinterface.h>
#include <assert.h>

using namespace std;

namespace grale
{

float calculateKappaThresholdFitness(const ProjectedImagesInterface &iface, const vector<int> &sourceIndices,
		                             const vector<float> &kappaThresholds)
{
	float fitness = 0;

	assert(sourceIndices.size() == kappaThresholds.size());

	for (size_t sIdx = 0 ; sIdx < sourceIndices.size() ; sIdx++)
	{
		const int s = sourceIndices[sIdx];
		const float threshold = kappaThresholds[sIdx];

		assert(s >= 0 && s < iface.getNumberOfSources());
		assert(iface.getNumberOfImages(s) == 1);
		assert(iface.getNumberOfImagePoints(s, 0) > 0);

		const int numPoints = iface.getNumberOfImagePoints(s, 0);
		const float *pKappas = iface.getConvergence(s, 0);

		for (int i = 0 ; i < numPoints ; i++)
		{
			if (pKappas[i] > threshold)
				fitness += (pKappas[i]-threshold);
		}
	}
	assert(!isnan(fitness));
	return fitness;
}

} // end namespace
