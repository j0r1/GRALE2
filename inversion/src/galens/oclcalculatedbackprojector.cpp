#include "oclcalculatedbackprojector.h"
#include "imagesdataextended.h"
#include "constants.h"
#include <memory>

namespace grale
{

using namespace std;
using namespace errut;

OclCalculatedBackProjector::OclCalculatedBackProjector()
{
}

OclCalculatedBackProjector::~OclCalculatedBackProjector()
{
}

bool_t OclCalculatedBackProjector::init(const std::vector<ImagesDataExtended *> &images, double angularScale)
{
	// Keep namespace clean
	storeOriginalData(images);

	if (images.size() == 0)
		return "No images data sets were specified";

	const int testSize = 10;
	Vector2Df arr1[testSize];
	float arr2[testSize*2];

	if (sizeof(arr1) != sizeof(arr2))
		return "Alignment error, can't use backprojector for OpenCL buffer";

	m_thetas.resize(images.size());
	m_angularScale = angularScale;

	m_sourceOffsets.clear();	
	int curOffset = 0;
	for (auto pImgDat : images)
	{
		m_sourceOffsets.push_back(curOffset);
		assert(pImgDat);
		int numImages = pImgDat->getNumberOfImages();
		for (int j = 0 ; j < numImages ; j++)
		{
			int numPoints = pImgDat->getNumberOfImagePoints(j);
			curOffset += numPoints;
		}
	}
		
	for (size_t s = 0 ; s < images.size() ; s++)
	{
		ImagesData *pImg = images[s];
		assert(pImg);

		int numImages = pImg->getNumberOfImages();

		assert(s < m_numTotalPoints.size());
		m_thetas[s].resize(m_numTotalPoints[s]);

		int pos = 0;
		for (int j = 0 ; j < numImages ; j++)
		{
			int numPoints = pImg->getNumberOfImagePoints(j);

			for (int k = 0 ; k < numPoints ; k++, pos++)
			{
				Vector2Dd theta = pImg->getImagePointPosition(j,k);
				theta /= m_angularScale;
				
				m_thetas[s][pos] = Vector2Df(theta.getX(), theta.getY());
			}			
		}
	}

	return true;
}

} // end namespace

