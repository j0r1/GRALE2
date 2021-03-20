#pragma once

#include "graleconfig.h"
#include "vector2d.h"
#include <errut/booltype.h>
#include <vector>

namespace grale
{

class ImagesDataExtended;
class ProjectedImagesInterface;

// ImagesData should have three images, on a grid, with same number of points:
// 0 - the central point
// 1 - offset in x direction
// 2 - offset in y direction
class NumericGradientCalculator
{
public:
	NumericGradientCalculator();
	~NumericGradientCalculator();
	errut::bool_t check(const ImagesDataExtended &img);
	Vector2Df getGradientFloat(int srcIdx, int ptIdx, const ProjectedImagesInterface &iface);
	Vector2Dd getGradientDouble(int srcIdx, int ptIdx, const ProjectedImagesInterface &iface);
	const std::vector<Vector2Df> &getGradientFloats(int srcIdx, const ProjectedImagesInterface &iface);
	const std::vector<Vector2Dd> &getGradientDoubles(int srcIdx, const ProjectedImagesInterface &iface);
	const std::vector<float> &getGradientSizes(int srcIdx, const ProjectedImagesInterface &iface);
	const std::vector<float> &getGradientSquaredSizes(int srcIdx, const ProjectedImagesInterface &iface);
private:
	std::vector<Vector2Df> m_tmpGradientsFloat;
	std::vector<Vector2Dd> m_tmpGradientsDouble;
	std::vector<float> m_tmpGradientsSizes;
};

} // namespace grale
