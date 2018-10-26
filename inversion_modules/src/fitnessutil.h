#ifndef GRALE_FITNESSUTIL_H

#define GRALE_FITNESSUTIL_H

#include <grale/triangleindices.h>
#include <grale/vector2d.h>
#include <vector>

namespace grale
{

class ProjectedImagesInterface;
class PointGroupStorage;
class FitnessComponentCache;

float getScaleFactor_PointImages(const ProjectedImagesInterface &interface,
		                         const std::vector<int> &sourceIndices,
								 const std::vector<float> &sourceDistanceFractions);

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const std::vector<int> &sourceIndices,
										  const std::vector<float> &sourceDistanceFractions,
										  float scale);

float calculateOverlapFitness_Extended(const PointGroupStorage &pointGroups, const ProjectedImagesInterface &iface,
		                               const std::vector<int> &sourceIndices,
									   const std::vector<bool> &rectFlags,
									   const std::vector<bool> &groupFlags,
									   FitnessComponentCache *pCache = 0);

float calculateNullFitness_PointImages(const ProjectedImagesInterface &interface, 
		                               const std::vector<int> &sourceIndices,
									   const std::vector<int> &nullIndices,
									   const std::vector<std::vector<TriangleIndices> > &nullTriangles,
									   const std::vector<float> &nullWeights);

float calculateNullFitness_ExtendedImages(const ProjectedImagesInterface &iface, 
		                               const std::vector<int> &sourceIndices,
									   const std::vector<int> &nullIndices,
									   const std::vector<std::vector<TriangleIndices> > &nullTriangles,
									   const std::vector<std::vector<double> > &nullTriangleAreas,
									   const std::vector<float> &nullWeights,
									   FitnessComponentCache *pCache = 0);

float calculateWeakLensingFitness(const ProjectedImagesInterface &interface, const std::vector<int> &weakIndices, 
							      const std::vector<bool> &reduced, const std::vector<double> &oneMinusKappaThreshold);

float calculateTimeDelayFitness(const ProjectedImagesInterface &iface, const std::vector<int> &sourceIndices);

float calculateTimeDelayFitnessExperimental(const ProjectedImagesInterface &iface, const std::vector<int> &sourceIndices);

float calculateKappaThresholdFitness(const ProjectedImagesInterface &iface, const std::vector<int> &sourceIndices,
		                             const std::vector<float> &kappaThresholds);

float calculateCausticPenaltyFitness(const ProjectedImagesInterface &iface,
		const std::vector<int> &sourceIndices,
		const std::vector<int> &gridIndices,
		const std::vector<std::vector<std::vector<std::pair<int,int> > > > &lineSegments,
		const std::vector<std::vector<std::vector<TriangleIndices> > > &critTriangles,
		std::vector<std::vector<std::vector<bool> > > &lineSegmentFlags,
		std::vector<std::vector<std::vector<Vector2D<float> > > > &lineSegmentIntersections,
		FitnessComponentCache *pCache);

} // end namespace

#endif // GRALE_FITNESSUTIL_H
