#include "fitnesscomponent_null.h"
#include "fitnessutil.h"
#include <grale/imagesdataextended.h>

#include <iostream>

using namespace std;

namespace grale
{

// FitnessComponent_NullSpacePointImages

FitnessComponent_NullSpacePointImages::FitnessComponent_NullSpacePointImages(FitnessComponentCache *pCache) 
	: FitnessComponent("pointimagenull", pCache)
{
	addRecognizedTypeName("pointimages");
	addRecognizedTypeName("pointnullgrid");
	addRecognizedTypeName("singlyimagedpoints");
	addRecognizedTypeName("extendedimages");
	addRecognizedTypeName("pointgroupimages");

	m_lastImgIdx = -1;
	m_lastImgDds = -1;
	m_lastImgDs = -1;
	m_lastImgNumImgs = 0;
}

FitnessComponent_NullSpacePointImages::~FitnessComponent_NullSpacePointImages()
{
}

bool FitnessComponent_NullSpacePointImages::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	const int numImg = imgDat.getNumberOfImages();

	auto finalize = [numImg, idx, this, &imgDat, &needCalcDeflections](const ImagesDataExtended &grpDat)
	{
		needCalcDeflections = true;

		m_lastImgIdx = idx;
		m_lastImgDds = imgDat.getDds();
		m_lastImgDs = imgDat.getDs();
		m_lastImgNumImgs = numImg;
		m_lastPointGroupImg = grpDat;
	};

	imgDat.getExtraParameter("type", typeName);
	if (typeName == "pointimages") 
	{
		if (numImg < 1)
		{
			setErrorString("Images data set doesn't contain any images");
			return false;
		}

		if (!isPointImage(imgDat))
		{
			setErrorString("Image is not a point image");
			return false;
		}
		
		string errStr;
		auto grpImg = addGroupsToPointImages(imgDat, errStr);
		if (!grpImg.get())
		{
			setErrorString("Couldn't add group info to point image: " + errStr);
			return false;
		}

		finalize(*grpImg.get());
		return true;
	}

	if (typeName == "singlyimagedpoints")
	{
		if (numImg != 1)
		{
			setErrorString("For the specified images data set, only a single 'image' is allowed");
			return false;
		}

		int numPoints = imgDat.getNumberOfImagePoints(0);
		if (numPoints < 1)
		{
			setErrorString("At least one point should be present in the 'image'");
			return false;
		}

		finalize(ImagesDataExtended()); // No groups are needed here, just add empty
		return true;
	}

	if (typeName == "extendedimages" || typeName == "pointgroupimages")
	{
		if (numImg < 1)
		{
			setErrorString("Images data set doesn't contain any images");
			return false;
		}

		if (imgDat.getNumberOfGroups() > 0) // ok, point groups are specified
		{
			finalize(imgDat);
			return true;
		}
		
		// no point groups, perhaps it's a point image and we can use this
		if (!isPointImage(imgDat))
		{
			cerr << "INFO: Image does not contain point groups and is not a point image" << endl;
			return true; // Just ignore this set
		}

		// Ok, it's a point image, so add group info
		string errStr;
		auto grpImg = addGroupsToPointImages(imgDat, errStr);
		if (!grpImg.get())
		{
			setErrorString("Couldn't add group info to point image: " + errStr);
			return false;
		}

		finalize(*grpImg.get());
		return true;
	}

	if (typeName != "pointnullgrid") // ignore the rest
		return true;

	if (m_lastImgIdx < 0)
	{
		setErrorString("No previous point image information that this null space grid can refer to");
		return false;
	}

	// TODO: use another check? Dds doesn't make much sense in a multi-plane scenario
	if (m_lastImgDds != imgDat.getDds() || m_lastImgDs != imgDat.getDs())
	{
		setErrorString("Null space grid distances are not the same as the ones from the last recorded point images");
		return false;
	}

	needCalcDeflections = true;

	double weight = 1.0f;

	if (imgDat.hasExtraParameter("weight") && !imgDat.getExtraParameter("weight", weight))
	{
		setErrorString("Null space grid specifies a weight, but seems to not be a real value: " + imgDat.getErrorString());
		return false;
	}

	if (weight < 0)
	{
		setErrorString("Specified weight for null space grid is negative");
		return false;
	}

	if (numImg != 1)
	{
		setErrorString("Null space data should contain only one image");
		return false;
	}

	if (!imgDat.hasTriangulation())
	{
		setErrorString("Null space data doesn't contain a triangulation");
		return false;
	}

	m_nullWeights.push_back((float)weight);
	
	int nullIdx = m_nullTriangles.size();
	m_nullTriangles.resize(nullIdx+1);
	m_nullTriangles[nullIdx].clear();

	vector<TriangleIndices> t;
	if (!imgDat.getTriangles(0, t))
	{
		setErrorString("Unable to obtain triangulation data from the null space data: " + imgDat.getErrorString());
		return false;
	}
	
	m_nullTriangles[nullIdx].resize(t.size());

	auto triangIt = t.begin();
	int i = 0;
	for ( ; triangIt != t.end() ; triangIt++, i++)
		m_nullTriangles[nullIdx][i] = *triangIt;

	m_sourceIndices.push_back(m_lastImgIdx);
	m_nullIndices.push_back(idx);
	m_pointGroups.add(m_lastPointGroupImg);

	addImagesDataIndex(m_lastImgIdx);
	addImagesDataIndex(idx);
	addToUsedImagesCount(m_lastImgNumImgs);
	//Note: disabled this again, so we can easilye two null spaces (perhaps with diffent resolution)
	//      for the same image
	//m_lastImgIdx = -1; // reset the referenced image

	return true;
}

bool FitnessComponent_NullSpacePointImages::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateNullFitness_PointImages(m_pointGroups, iface, m_sourceIndices, 
			                                   m_nullIndices, m_nullTriangles, m_nullWeights);
	return true;
}

// FitnessComponent_NullSpaceExtendedImages

FitnessComponent_NullSpaceExtendedImages::FitnessComponent_NullSpaceExtendedImages(FitnessComponentCache *pCache) 
	: FitnessComponent("extendedimagenull", pCache)
{
	addRecognizedTypeName("extendedimages");
	addRecognizedTypeName("extendednullgrid");

	m_lastImgIdx = -1;
	m_lastImgDds = -1;
	m_lastImgDs = -1;
	m_lastImgNumImgs = 0;
}

FitnessComponent_NullSpaceExtendedImages::~FitnessComponent_NullSpaceExtendedImages()
{
}

bool FitnessComponent_NullSpaceExtendedImages::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName == "extendedimages") 
	{
		int numImg = imgDat.getNumberOfImages();
		if (numImg < 1)
		{
			setErrorString("Images data set doesn't contain any images");
			return false;
		}

		if (numImg == 1)
		{
			if (imgDat.getNumberOfImagePoints(0) != 1)
			{
				setErrorString("In case only a single image is present, it must be a single point (used for null space)");
				return false;
			}
		}
		else
		{
			for (int i = 0 ; i < numImg ; i++)
			{
				int numPoints = imgDat.getNumberOfImagePoints(i);

				if (numPoints < 1)
				{
					setErrorString("An empty image was present");
					return false;
				}
			}
		}
		
		// Just record some information at this point

		needCalcDeflections = true;

		m_lastImgIdx = idx;
		m_lastImgDds = imgDat.getDds();
		m_lastImgDs = imgDat.getDs();
		m_lastImgNumImgs = numImg;

		return true;
	}

	if (typeName != "extendednullgrid") // ignore the rest
		return true;

	if (m_lastImgIdx < 0)
	{
		setErrorString("No previous extended image information that this null space grid can refer to");
		return false;
	}

	// TODO: use another check? Dds doesn't make much sense in a multi-plane scenario
	if (m_lastImgDds != imgDat.getDds() || m_lastImgDs != imgDat.getDs())
	{
		setErrorString("Null space grid distances are not the same as the ones from the last recorded extended images");
		return false;
	}

	needCalcDeflections = true;

	double weight = 1.0f;

	if (imgDat.hasExtraParameter("weight") && !imgDat.getExtraParameter("weight", weight))
	{
		setErrorString("Null space grid specifies a weight, but seems to not be a real value: " + imgDat.getErrorString());
		return false;
	}

	if (weight < 0)
	{
		setErrorString("Specified weight for null space grid is negative");
		return false;
	}

	int numImg = imgDat.getNumberOfImages();
	if (numImg != 1)
	{
		setErrorString("Null space data should contain only one image");
		return false;
	}

	if (!imgDat.hasTriangulation())
	{
		setErrorString("Null space data doesn't contain a triangulation");
		return false;
	}

	m_nullWeights.push_back((float)weight);
	
	int nullIdx = m_nullTriangles.size();
	m_nullTriangles.resize(nullIdx+1);
	m_nullTriangleAreas.resize(nullIdx+1);

	vector<TriangleIndices> t;
	if (!imgDat.getTriangles(0, t))
	{
		setErrorString("Unable to obtain triangulation data from the null space data: " + imgDat.getErrorString());
		return false;
	}
	
	m_nullTriangles[nullIdx].resize(t.size());
	m_nullTriangleAreas[nullIdx].resize(t.size());

	auto triangIt = t.begin();
	size_t i = 0;
	for ( ; triangIt != t.end() ; triangIt++, i++)
		m_nullTriangles[nullIdx][i] = *triangIt;

	vector<Vector2D<double> > points;
	points.resize(imgDat.getNumberOfImagePoints(0));

	for (i = 0 ; i < points.size() ; i++)
		points[i] = imgDat.getImagePointPosition(0, i);

	for (i = 0 ; i < m_nullTriangles[nullIdx].size() ; i++)
	{
		Vector2D<double> p1 = points[m_nullTriangles[nullIdx][i].getIndex(0)];
		Vector2D<double> p2 = points[m_nullTriangles[nullIdx][i].getIndex(1)];
		Vector2D<double> p3 = points[m_nullTriangles[nullIdx][i].getIndex(2)];

		Triangle2D<double> t(p1,p2,p3);
		m_nullTriangleAreas[nullIdx][i] = t.getArea()/(ANGLE_ARCSEC*ANGLE_ARCSEC);
	}

	m_sourceIndices.push_back(m_lastImgIdx);
	m_nullIndices.push_back(idx);

	addImagesDataIndex(m_lastImgIdx);
	addImagesDataIndex(idx);
	addToUsedImagesCount(m_lastImgNumImgs);
	//Note: disabled this again, so we can easilye two null spaces (perhaps with diffent resolution)
	//      for the same image
	//m_lastImgIdx = -1; // reset the referenced image

	return true;
}

bool FitnessComponent_NullSpaceExtendedImages::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateNullFitness_ExtendedImages(iface, m_sourceIndices, m_nullIndices, m_nullTriangles,
			                                      m_nullTriangleAreas, m_nullWeights);
	return true;
}

} // end namespace
