#include "numericgradientcalculator.h"
#include <grale/imagesdataextended.h>
#include <grale/projectedimagesinterface.h>
#include <sstream>

using namespace std;
using namespace errut;

namespace grale
{

NumericGradientCalculator::NumericGradientCalculator()
{

}

NumericGradientCalculator::~NumericGradientCalculator()
{

}

bool_t NumericGradientCalculator::check(const ImagesDataExtended &img)
{
	int numImg = img.getNumberOfImages();

	if (numImg != 3)
		return "Three images are required in the data set (" + to_string(numImg) + " are present), one for the base points, one for the dx points, and one for the dy points";
	
	int numPts = img.getNumberOfImagePoints(0);
	if (numPts == 0)
		return "No points are present for the base image (index 0)";

	int numPts1 = img.getNumberOfImagePoints(1);
	int numPts2 = img.getNumberOfImagePoints(2);

	if (numPts1 != numPts)
		return "The dx image (index 1) does not contain " + to_string(numPts) + " points, but " + to_string(numPts1);
	if (numPts2 != numPts)
		return "The dy image (index 2) does not contain " + to_string(numPts) + " points, but " + to_string(numPts2);

	auto to_string_precise = [](double x)
	{
		stringstream ss;
		ss.precision(10);
		ss << x;
		return ss.str();
	};

	auto checkDirection = [&img,numPts,to_string_precise](int otherImgIdx) -> bool_t
	{
		assert(otherImgIdx == 1 || otherImgIdx == 2);
		int thisComp = (otherImgIdx == 1)?0:1; // 0 for x, 1 for y
		int otherComp = 1-thisComp; // switches them
		string componentName[] { "X", "Y" };
		string imgName[] { "base", "dx", "dy" };

		for (int i = 0 ; i < numPts ; i++)
		{
			Vector2Dd p0 = img.getImagePointPosition(0, i);
			Vector2Dd p1 = img.getImagePointPosition(otherImgIdx, i);
			double diffOther = p1.getComponents()[otherComp] - p0.getComponents()[otherComp];
			double diffThis = p1.getComponents()[thisComp] - p0.getComponents()[thisComp];
			if (diffOther != 0)
				return "The " + componentName[otherComp] + " value for the " + imgName[otherImgIdx] +
				       " image is not the same as for the base one (point " + to_string(i) + ") (diff is " + to_string_precise(diffOther) + ")";

			if (diffThis <= 0)
				return "The " + componentName[thisComp] + " value for the " + imgName[otherImgIdx] +
				       " image should be larger than zero, but is " + to_string_precise(diffThis) + " (point " + to_string(i) + ")";
		}

		return true;
	};

	bool_t r;
	if (!(r = checkDirection(1)))
		return r;
	if (!(r = checkDirection(2)))
		return r;

	// TODO: also check the size of the differences? How large/small should they be?
	return true;
}

Vector2Dd NumericGradientCalculator::getGradientDouble(int srcIdx, int ptIdx, const ProjectedImagesInterface &iface)
{
	assert(srcIdx >= 0 && srcIdx < iface.getNumberOfSources());
	assert(iface.getNumberOfImages(srcIdx) == 3);
	assert(iface.getNumberOfImagePoints(srcIdx, 0) == iface.getNumberOfImagePoints(srcIdx, 1));
	assert(iface.getNumberOfImagePoints(srcIdx, 0) == iface.getNumberOfImagePoints(srcIdx, 2));
	assert(iface.getNumberOfImagePoints(srcIdx, 0) > 0);
	assert(ptIdx >= 0 && ptIdx < iface.getNumberOfImagePoints(srcIdx, 0));

	float v0 = iface.getConvergence(srcIdx, 0)[ptIdx];
	float v1 = iface.getConvergence(srcIdx, 1)[ptIdx];
	float v2 = iface.getConvergence(srcIdx, 2)[ptIdx];
	Vector2Df x0 = iface.getThetas(srcIdx, 0)[ptIdx];
	Vector2Df x1 = iface.getThetas(srcIdx, 1)[ptIdx];
	Vector2Df x2 = iface.getThetas(srcIdx, 2)[ptIdx];
	double dx = x1.getX()-x0.getX();
	double dy = x2.getY()-x0.getY();
	double dVx = v1-v0;
	double dVy = v2-v0;
	return Vector2Dd(dVx/dx, dVy/dy);
}

Vector2Df NumericGradientCalculator::getGradientFloat(int srcIdx, int ptIdx, const ProjectedImagesInterface &iface)
{
	Vector2Dd g = getGradientDouble(srcIdx, ptIdx, iface);
	return Vector2Df((float)g.getX(), (float)g.getY());
}

template<class T, typename F>
const vector<T> &fillInVector(vector<T> &v, int srcIdx, const ProjectedImagesInterface &iface, F getFunction)
{
	assert(srcIdx >= 0 && srcIdx < iface.getNumberOfSources());
	assert(iface.getNumberOfImages(srcIdx) == 3);

	const int numPoints = iface.getNumberOfImagePoints(srcIdx, 0);
	v.resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
	{
		// TODO: optimize this, shouldn't be a need to calc dx/dy each time
		v[i] = getFunction(srcIdx, i, iface);
	}
	return v;
}

const vector<Vector2Dd> &NumericGradientCalculator::getGradientDoubles(int srcIdx, const ProjectedImagesInterface &iface)
{
	return fillInVector(m_tmpGradientsDouble, srcIdx, iface,
	[this](int srcIdx, int ptIdx, auto &iface) {
		return getGradientDouble(srcIdx, ptIdx, iface);
	});
}

const vector<Vector2Df> &NumericGradientCalculator::getGradientFloats(int srcIdx, const ProjectedImagesInterface &iface)
{
	return fillInVector(m_tmpGradientsFloat, srcIdx, iface,
	[this](int srcIdx, int ptIdx, auto &iface) {
		return getGradientFloat(srcIdx, ptIdx, iface);
	});
}

const std::vector<float> &NumericGradientCalculator::getGradientSizes(int srcIdx, const ProjectedImagesInterface &iface)
{
	assert(srcIdx >= 0 && srcIdx < iface.getNumberOfSources());
	assert(iface.getNumberOfImages(srcIdx) == 3);

	const vector<Vector2Dd> &v = getGradientDoubles(srcIdx, iface);
	m_tmpGradientsSizes.resize(v.size());
	
	for (size_t i = 0 ; i < m_tmpGradientsSizes.size() ; i++)
		m_tmpGradientsSizes[i] = (float)v[i].getLength();

	return m_tmpGradientsSizes;
}

}
