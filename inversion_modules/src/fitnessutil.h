#ifndef GRALE_FITNESSUTIL_H

#define GRALE_FITNESSUTIL_H

#include <grale/triangleindices.h>
#include <grale/vector2d.h>
#include <grale/imagesdataextended.h>
#include <vector>
#include <memory>

namespace grale
{

class ProjectedImagesInterface;
class PointGroupStorage;
class FitnessComponentCache;
class NumericGradientCalculator;

std::shared_ptr<ImagesDataExtended> addGroupsToPointImages(const ImagesDataExtended &imgDat, std::string &errStr);

float getScaleFactor_PointImages(const ProjectedImagesInterface &interface,
		                         const std::vector<int> &sourceIndices);

class ScaleFactorWorkspace
{
public:
	std::vector<float> m_floats;
	std::vector<std::vector<float>> m_vecFloat;
};

enum PointImageScaleType { MinMax, MAD };
void getScaleFactors_PointImages(const ProjectedImagesInterface &interface,
		                         const std::vector<int> &sourceIndices,
								 const std::vector<int> &sourceGroup,
								 std::vector<float> &scaleFactors, // one for each group
								 ScaleFactorWorkspace &ws,
								 PointImageScaleType scaleType);

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const std::vector<int> &sourceIndices,
										  float scale);

float calculateOverlapFitness_PointImages(const ProjectedImagesInterface &interface, 
		                                  const std::vector<int> &sourceIndices,
										  const std::vector<int> &sourceGroups,
										  const std::vector<float> &scaleFactors);

enum PointGroupRMSType { AllBetas, AverageBeta };
float calculateOverlapFitness_PointGroups(const PointGroupStorage &pointGroups,
		                                  const ProjectedImagesInterface &interface,
										  const std::vector<int> &sourceIndices,
										  PointGroupRMSType t);

float calculateOverlapFitness_Extended(const PointGroupStorage &pointGroups, const ProjectedImagesInterface &iface,
		                               const std::vector<int> &sourceIndices,
									   const std::vector<bool> &rectFlags,
									   const std::vector<bool> &groupFlags,
									   FitnessComponentCache *pCache = 0);

float calculateNullFitness_PointImages(const PointGroupStorage &pointGroups,
		                               const ProjectedImagesInterface &interface, 
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

enum WeakLensingType { RealShear, RealReducedShear, AveragedEllipticities, BayesianEllipticities };
float calculateWeakLensingFitness(const ProjectedImagesInterface &interface, const std::vector<int> &weakIndices, 
							      WeakLensingType weakType, const std::vector<float> &oneMinusKappaThreshold);

float calculateTimeDelayFitnessPaper2009(const ProjectedImagesInterface &iface, const std::vector<int> &sourceIndices,
		                        const std::vector<float> &tdScaleFactors = std::vector<float>());

float calculateTimeDelayFitnessNoSrc(const ProjectedImagesInterface &iface, const std::vector<int> &sourceIndices,
                                             const std::vector<float> &tdScaleFactors = std::vector<float>());

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
