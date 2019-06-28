#include "fitnesscomponent.h"
#include "fitnessutil.h"
#include <grale/imagesdataextended.h>
#include <grale/projectedimagesinterface.h>
#include <limits>
#include <list>
#include <sstream>

#include <iostream>

using namespace std;

namespace grale
{

FitnessComponentCache::FitnessComponentCache(int num)
{
	m_numImgDatas = num;
	m_estimatedSourceShapes.resize(num);
	m_sourceScales.resize(num);
	m_clear = true;
}

FitnessComponentCache::~FitnessComponentCache()
{
	clear();
}

void FitnessComponentCache::clear()
{
	if (m_clear)
		return;

	for (int i = 0 ; i < m_numImgDatas ; i++)
	{
		m_estimatedSourceShapes[i].m_isEstSrcSet = false;
		m_sourceScales[i].m_isSet = false;
	}

	m_clear = true;
}

void FitnessComponentCache::setEstimatedSource(int idx, const Polygon2D<float> &poly, 
                                               float minX, float maxX, float minY, float maxY)
{
	assert(idx >= 0 && idx < m_estimatedSourceShapes.size());
	SourceShape &s = m_estimatedSourceShapes[idx];

	s.m_isEstSrcSet = true;
	s.m_estimatedSource = poly;
	s.m_srcMinX = minX;
	s.m_srcMaxX = maxX;
	s.m_srcMinY = minY;
	s.m_srcMaxY = maxY;

	m_clear = false;
}

bool FitnessComponentCache::getEstimatedSource(int idx, Polygon2D<float> &poly, 
                                               float &minX, float &maxX, float &minY, float &maxY)
{
	assert(idx >= 0 && idx < m_estimatedSourceShapes.size());
	SourceShape &s = m_estimatedSourceShapes[idx];

	if (!s.m_isEstSrcSet)
		return false;

	poly = s.m_estimatedSource;
	minX = s.m_srcMinX;
	maxX = s.m_srcMaxX;
	minY = s.m_srcMinY;
	maxY = s.m_srcMaxY;

	return true;
}

void FitnessComponentCache::setEstimatedSourceScale(int idx, float scale)
{
	assert(idx >= 0 && idx < m_sourceScales.size());

	m_sourceScales[idx].m_isSet = true;
	m_sourceScales[idx].m_scale = scale;
	m_clear = false;
}

bool FitnessComponentCache::getEstimatedSourceScale(int idx, float &scale)
{
	assert(idx >= 0 && idx < m_sourceScales.size());
	if (!m_sourceScales[idx].m_isSet)
		return false;

	scale = m_sourceScales[idx].m_scale;
	return true;
}

FitnessComponent::FitnessComponent(const string &name, FitnessComponentCache *pCache) : errut::ErrorBase(name)
{
	m_pCache = pCache;
	m_priority = numeric_limits<int>::max();
	m_priorityIndex = -1;
	m_usedImagesCount = 0;
}

FitnessComponent::~FitnessComponent()
{
}

bool FitnessComponent::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	setErrorString("FitnessComponent::inspectImagesData needs to be implemented in a subclass");
	return false;
}

bool FitnessComponent::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	setErrorString("FitnessComponent::calculateFitness needs to be implemented in a subclass");
	return false;
}

void FitnessComponent::addImagesDataIndex(int idx)
{
	assert(idx >= 0);
	m_usedImagesDataIndices.push_back(idx);
}

void FitnessComponent::addRecognizedTypeName(const string &n)
{
	m_recognizedTypeNames.push_back(n);
}

void FitnessComponent::addToUsedImagesCount(int c)
{
	assert(c >= 0);
	m_usedImagesCount += c;
}

bool FitnessComponent::isPointImage(const ImagesData &imgDat)
{
	for (int i = 0 ; i < imgDat.getNumberOfImages() ; i++)
	{
		if (imgDat.getNumberOfImagePoints(i) != 1)
			return false;
	}
	return true;
}

// FitnessComponent_PointImagesOverlap

FitnessComponent_PointImagesOverlap::FitnessComponent_PointImagesOverlap(FitnessComponentCache *pCache) 
	: FitnessComponent("pointimageoverlap", pCache)
{
	addRecognizedTypeName("pointimages");
	m_scaleType = MinMax;
}

FitnessComponent_PointImagesOverlap::~FitnessComponent_PointImagesOverlap()
{
}

bool FitnessComponent_PointImagesOverlap::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
							   bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
							   bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
							   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "pointimages") // ignore
		return true;
	
	int numImg = imgDat.getNumberOfImages();
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

	double Dds = imgDat.getDds();
	double Ds = imgDat.getDs();

	if (Dds <= 0 || Ds <= 0)
	{
		setErrorString("Source/lens and source/observer distances must be positive");
		return false;
	}

	// If there's only a single point, it's not a constraint for this fitness component
	// (but can still be used as a null space constraint)
	if (numImg == 1)
		return true;

	needCalcDeflections = true;

	float distFrac = (float)(Dds/Ds);
	m_distanceFractions.push_back(distFrac);

	addImagesDataIndex(idx);
	addToUsedImagesCount(numImg);

	if (imgDat.hasExtraParameter("groupname"))
	{
		string v;
		if (!imgDat.getExtraParameter("groupname", v))
		{
			setErrorString("'groupname' parameter is present, but is not a string");
			return false;
		}

		if (v.length() == 0)
		{
			setErrorString("'groupname' parameter is present, but the value is an empty string");
			return false;
		}

		m_groupnameIndices[v].push_back(idx);

		// TODO: check for same distance fraction?
	}
	else
		m_nogroupnameIndices.push_back(idx);

	if (imgDat.hasExtraParameter("usescalefrom"))
	{
		string v;
		if (!imgDat.getExtraParameter("usescalefrom", v))
		{
			setErrorString("'usescalefrom' is present, but is not a string");
			return false;
		}

		if (v.length() == 0)
		{
			setErrorString("'usescalefrom' is present, but the value is an empty string");
			return false;
		}

		m_useScaleNames[idx] = v;
	}

	return true;
}

bool FitnessComponent_PointImagesOverlap::finalize()
{
	if (m_setScaleGroups.size() == m_distanceFractions.size())
		return true;

	//cerr << "DEBUG: finalize in FitnessComponent_PointImagesOverlap, for " << (void*)this << endl;

	// Need to build the group indices from the map
	m_setScaleGroups.resize(m_distanceFractions.size(), -1);

	const auto &imageIndices = getUsedImagesDataIndices();
	map<int,int> imgIndexToArrayIndex;
	for (int i = 0 ; i < imageIndices.size() ; i++)
		imgIndexToArrayIndex[imageIndices[i]] = i;

	map<string, int> groupNameToIndex;
	int nextGroup = 0;

	auto addIndices = [&groupNameToIndex,&nextGroup,&imgIndexToArrayIndex,this](const vector<int> &indices, const string &name)
	{
		for (auto idx : indices)
		{
			assert(imgIndexToArrayIndex.find(idx) != imgIndexToArrayIndex.end());
			int i = imgIndexToArrayIndex[idx];

			assert(i >= 0 && i < m_setScaleGroups.size());
			m_setScaleGroups[i] = nextGroup;
		}
		groupNameToIndex[name] = nextGroup;
	};

	// Unnamed points go into group zero
	addIndices(m_nogroupnameIndices, "");
	nextGroup++;

	for (auto it : m_groupnameIndices)
	{
		addIndices(it.second, it.first);
		nextGroup++;
	}

	// check that every point has a group number (was initialized to -1)
	for (auto idx : m_setScaleGroups) { assert(m_setScaleGroups[idx] >= 0); }

	// Reserve space in scale factor array
	m_scaleFactors.resize(nextGroup);

	// To allow a point to use a scale from another group, the m_useScaleGroups
	// array is used
	// Initialize this to m_setScaleGroups, and override when specified
	m_useScaleGroups = m_setScaleGroups;

	for (auto it : m_useScaleNames)
	{
		int imageIndex = it.first;
		string groupName = it.second;

		auto it2 = groupNameToIndex.find(groupName);
		if (it2 == groupNameToIndex.end())
		{
			setErrorString("Group name '" + groupName + "' was mentioned in a 'usescalefrom' setting, but is not present as a 'groupname' setting");
			return false;
		}

		int groupIndex = it2->second;

		assert(imageIndex >= 0 && imageIndex < m_useScaleGroups.size());
		assert(groupIndex >= 0 && groupIndex < m_scaleFactors.size());
		m_useScaleGroups[imageIndex] = groupIndex;
	}

	//cerr << "m_setScaleGroups: " << m_setScaleGroups.size() << endl;
	//cerr << "m_useScaleGroups: " << m_useScaleGroups.size() << endl;
	//cerr << "m_distanceFractions: " << m_distanceFractions.size() << endl;
	return true;
}

bool FitnessComponent_PointImagesOverlap::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	//cerr << "In calculateFitness " << (void*)this << endl;
	//cerr << "m_setScaleGroups: " << m_setScaleGroups.size() << endl;
	//cerr << "m_useScaleGroups: " << m_useScaleGroups.size() << endl;
	//cerr << "m_distanceFractions: " << m_distanceFractions.size() << endl;
	assert(m_setScaleGroups.size() == m_distanceFractions.size());
	assert(m_useScaleGroups.size() == m_distanceFractions.size());

	getScaleFactors_PointImages(iface, getUsedImagesDataIndices(), m_distanceFractions, m_setScaleGroups, m_scaleFactors, m_workspace, m_scaleType);
	fitness = calculateOverlapFitness_PointImages(iface, getUsedImagesDataIndices(), m_distanceFractions, m_useScaleGroups, m_scaleFactors);
	return true;
}

// FitnessComponent_PointGroupOverlap

FitnessComponent_PointGroupOverlap::FitnessComponent_PointGroupOverlap(FitnessComponentCache *pCache) 
	: FitnessComponent("pointgroupoverlap", pCache)
{
	addRecognizedTypeName("pointgroupimages");
	m_rmsType = AllBetas;
}

FitnessComponent_PointGroupOverlap::~FitnessComponent_PointGroupOverlap()
{
}

bool FitnessComponent_PointGroupOverlap::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "pointgroupimages") // ignore
		return true;
	
	int numImg = imgDat.getNumberOfImages();
	if (numImg < 1)
	{
		setErrorString("Images data set doesn't contain any images");
		return false;
	}
	
	double Dds = imgDat.getDds();
	double Ds = imgDat.getDs();
	if (Dds <= 0 || Ds <= 0)
	{
		setErrorString("Source/lens and source/observer distances must be positive");
		return false;
	}
	float distFrac = (float)(Dds/Ds);
	m_distanceFractions.push_back(distFrac);

	if (imgDat.getNumberOfGroups() > 0)
		m_pointGroups.add(&imgDat);
	else // no point groups, is this a point image?
	{
		if (!isPointImage(imgDat))
		{
			setErrorString("No point groups present, and not a point image");
			return false;
		}

		string errStr;
		auto grpImg = addGroupsToPointImages(imgDat, errStr);
		if (!grpImg.get())
		{
			setErrorString("Couldn't add group info to point image: " + errStr);
			return false;
		}
		m_pointGroups.add(grpImg.get());
	}
	
	// Ok everything checks out

	needCalcDeflections = true;
	needCalcDeflDeriv = true;

	addImagesDataIndex(idx);
	addToUsedImagesCount(numImg);

	return true;
}

bool FitnessComponent_PointGroupOverlap::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateOverlapFitness_PointGroups(m_pointGroups, iface, getUsedImagesDataIndices(), 
			                                   m_distanceFractions, m_rmsType);
	return true;
}

// FitnessComponent_ExtendedImagesOverlap

FitnessComponent_ExtendedImagesOverlap::FitnessComponent_ExtendedImagesOverlap(FitnessComponentCache *pCache) 
	: FitnessComponent("extendedimageoverlap", pCache)
{
	addRecognizedTypeName("extendedimages");
	m_extendedSourceCount = 0;
}

FitnessComponent_ExtendedImagesOverlap::~FitnessComponent_ExtendedImagesOverlap()
{
}

bool FitnessComponent_ExtendedImagesOverlap::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "extendedimages") // ignore
		return true;
	
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
		int countPointImages = 0;
		for (int i = 0 ; i < numImg ; i++)
		{
			int numPoints = imgDat.getNumberOfImagePoints(i);
			if (numPoints == 0)
			{
				setErrorString("An empty image is present");
				return false;
			}
			else if (numPoints == 1)
				countPointImages++;
		}

		if (countPointImages != numImg)
			m_extendedSourceCount++;
	}

	double Dds = imgDat.getDds();
	double Ds = imgDat.getDs();

	if (Dds <= 0 || Ds <= 0)
	{
		setErrorString("Source/lens and source/observer distances must be positive");
		return false;
	}

	bool useRects = true;
	bool usePointGroups = (imgDat.getNumberOfGroups() > 0)?true:false;

	if (imgDat.hasExtraParameter("userectangles"))
	{
		if (!imgDat.getExtraParameter("userectangles", useRects))
		{
			setErrorString("Extra parameter 'userectangles' is present, but doesn't appear to be a boolean: " + imgDat.getErrorString());
			return false;
		}
	}

	if (imgDat.hasExtraParameter("usepointgroups"))
	{
		if (!imgDat.getExtraParameter("usepointgroups", usePointGroups))
		{
			setErrorString("Extra parameter 'usepointgroups' is present, but doesn't appear to be a boolean: " + imgDat.getErrorString());
			return false;
		}
	}

	if (usePointGroups && imgDat.getNumberOfGroups() == 0)
	{
		setErrorString("Point group overlap was requested, but images data set does not have any groups");
		return false;
	}

	if (!useRects && !usePointGroups) // TODO: should we just ignore this?
	{
		setErrorString("Neither surrounding rectangles nor point groups are used, image appears to be useless");
		return false;
	}

	// Ok everything checks out

	// If there's only a single point, it's not a constraint for this fitness component
	// (but can still be used as a null space constraint)
	if (numImg > 1)
	{
		needCalcDeflections = true;

		m_pointGroups.add(&imgDat);
		m_rectFlags.push_back(useRects);
		m_useGroupsFlags.push_back(usePointGroups);

		addImagesDataIndex(idx);
		addToUsedImagesCount(numImg);
	}

	return true;
}

bool FitnessComponent_ExtendedImagesOverlap::finalize()
{
	if (m_extendedSourceCount == 0)
	{
		setErrorString("At least one set of extended images needs to be present to use this fitness measure");
		return false;
	}
	return true;
}

bool FitnessComponent_ExtendedImagesOverlap::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateOverlapFitness_Extended(m_pointGroups, iface, getUsedImagesDataIndices(), 
			                                   m_rectFlags, m_useGroupsFlags, getCache());
	return true;
}

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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
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

		m_pointGroups.add(grpDat);
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
			setErrorString("Image does not contain point groups and is not a point image");
			return false;
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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
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

// FitnessComponent_WeakLensing

FitnessComponent_WeakLensing::FitnessComponent_WeakLensing(FitnessComponentCache *pCache) 
	: FitnessComponent("weaklensing", pCache)
{
	addRecognizedTypeName("sheardata");
	m_wlType = AveragedEllipticities;
}

FitnessComponent_WeakLensing::~FitnessComponent_WeakLensing()
{
}

bool FitnessComponent_WeakLensing::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "sheardata")
		return true; // ignore

	if (imgDat.getNumberOfImages() != 1)
	{
		setErrorString("Each images data instance can only contain one 'image' for shear info");
		return false;
	}

	if (!imgDat.hasShearInfo())
	{
		setErrorString("No shear info is present");
		return false;
	}

	needCalcDeflDeriv = true;
	needCalcShear = true;
	needCalcConvergence = true;
	storeOrigShear = true;

	double threshold = 0;
	if (!imgDat.getExtraParameter("threshold", threshold))
	{
		setErrorString("Shear data instance does not contain a (real valued) 'threshold' parameter for |1-kappa|");
		return false;
	}
	if (threshold < 0)
	{
		setErrorString("The threshold value for |1-kappa| must be positive or zero");
		return false;
	}

	m_thresholds.push_back((float)threshold);

	addImagesDataIndex(idx);

	return true;
}

bool FitnessComponent_WeakLensing::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateWeakLensingFitness(iface, getUsedImagesDataIndices(), m_wlType, m_thresholds);
	return true;
}

// FitnessComponent_TimeDelay

FitnessComponent_TimeDelay::FitnessComponent_TimeDelay(FitnessComponentCache *pCache) 
	: FitnessComponent("timedelay", pCache),
	  m_fitnessType(Paper2009),
	  m_relative(false)
{
	addRecognizedTypeName("pointimages");
	addRecognizedTypeName("extendedimages");
}

FitnessComponent_TimeDelay::~FitnessComponent_TimeDelay()
{
}

bool FitnessComponent_TimeDelay::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "extendedimages" && typeName != "pointimages")
		return true; // not relevant, ignore

	if (!imgDat.hasTimeDelays())
	{
		bool useTd = false;
		if (imgDat.hasExtraParameter("timedelay"))
		{
			if (!imgDat.getExtraParameter("timedelay", useTd))
			{
				setErrorString("Extra parameter 'timedelay' was specified, but is not a boolean");
				return false;
			}
			if (useTd)
			{
				setErrorString("Extra parameter 'timedelay' is set, but no time delays are present");
				return false;
			}
		}
		return true; // no time delays present, and usage was not requested. ignore
	}
	
	// Time delays are present in image, by default they should be used
	if (imgDat.hasExtraParameter("timedelay"))
	{
		bool useTd;

		if (!imgDat.getExtraParameter("timedelay", useTd))
		{
			setErrorString("Extra parameter 'timedelay' was specified, but is not a boolean");
			return false;
		}
		if (!useTd) // time delay info was explicitly specified to be skipped
			return true; // just ignore this
	}

	// Check that there are no negative time delay values: in the past this
	// was used to ignore a specific value, we'll generate an error if such
	// a value is still present
	int numTimeDelays = imgDat.getNumberOfTimeDelays();
	for (int tIdx = 0 ; tIdx < numTimeDelays ; tIdx++)
	{
		int img, point;
		double delay;
		imgDat.getTimeDelay(tIdx, &img, &point, &delay);
		if (delay < 0)
		{
			stringstream ss;
			ss << "A negative time delay (" << delay << ") was present for image " << img << ", point " << point;
			setErrorString(ss.str());
			return false;
		}
	}

	float tdScaleFactor = 0;
	if (imgDat.hasExtraParameter("timedelay_scalefactor"))
	{
		if (!(m_fitnessType == Paper2009 && m_relative == false))
		{
			setErrorString("'timedelay_scalefactor' was set but is only supported for the Paper2009 method and non-relative version");
			return false;
		}

		double v;

		if (!imgDat.getExtraParameter("timedelay_scalefactor", v))
		{
			setErrorString("Extra parameter 'timedelay_scalefactor' was specified, but is not a floating point value");
			return false;
		}
		if (v <= 0)
		{
			setErrorString("Extra parameter 'timedelay_scalefactor' should be positive");
			return false;
		}
		tdScaleFactor = (float)v;
	}
	
	needCalcPotential = true;
	storeOrigTimeDelay = true;

	//cerr << "Added TD fitness for " << idx << endl;
	addImagesDataIndex(idx);
	m_tdScaleFactors.push_back(tdScaleFactor);

	return true;
}

bool FitnessComponent_TimeDelay::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	if (m_relative && m_fitnessType != Paper2009)
	{
		setErrorString("Currently only the Paper2009 algorithm is supported for relative time delays");
		return false;
	}

	if (m_fitnessType == Paper2009)
	{
		if (!m_relative)
			fitness = calculateTimeDelayFitness(iface, getUsedImagesDataIndices(), m_tdScaleFactors);
		else
			fitness = calculateTimeDelayFitness_Relative(iface, getUsedImagesDataIndices(), m_referencePoints);
	}
	else if (m_fitnessType == ExpI)
		fitness = calculateTimeDelayFitnessExperimental(iface, getUsedImagesDataIndices());
	else if (m_fitnessType == ExpII)
		fitness = calculateTimeDelayFitnessExperimental2(iface, getUsedImagesDataIndices());
	else
	{
		setErrorString("Unknown TD fitness type");
		return false;
	}
	return true;
}

// FitnessComponent_KappaThreshold

FitnessComponent_KappaThreshold::FitnessComponent_KappaThreshold(FitnessComponentCache *pCache) 
	: FitnessComponent("kappathreshold", pCache)
{
	addRecognizedTypeName("kappathresholdpoints");
}

FitnessComponent_KappaThreshold::~FitnessComponent_KappaThreshold()
{
}

bool FitnessComponent_KappaThreshold::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "kappathresholdpoints")
		return true; // not relevant, ignore

	string thresStr;
	double threshold = 1.0;

	if (imgDat.hasExtraParameter("threshold"))
	{
		if (!imgDat.getExtraParameter("threshold", threshold))
		{
			setErrorString("Extra parameter 'threshold' is present, but doesn't seem to be a real value");
			return false;
		}
	}

	if (threshold < 0)
	{
		setErrorString("The specified threshold is negative");
		return false;
	}

	if (imgDat.getNumberOfImages() != 1)
	{
		setErrorString("Only a single large 'image' (containing a set of points for which kappa should be inspected) is allowed");
		return false;
	}
	if (imgDat.getNumberOfImagePoints(0) < 1)
	{
		setErrorString("No image points are present");
		return false;
	}

	needCalcDeflDeriv = true;
	needCalcConvergence = true;

	addImagesDataIndex(idx);
	m_thresholds.push_back((float)threshold);
	return true;
}

bool FitnessComponent_KappaThreshold::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateKappaThresholdFitness(iface, getUsedImagesDataIndices(), m_thresholds);
	return true;
}

// FitnessComponent_CausticPenalty

FitnessComponent_CausticPenalty::FitnessComponent_CausticPenalty(FitnessComponentCache *pCache)
	: FitnessComponent("causticpenalty", pCache)
{
	addRecognizedTypeName("causticgrid");
}

FitnessComponent_CausticPenalty::~FitnessComponent_CausticPenalty()
{
}

bool FitnessComponent_CausticPenalty::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName == "extendedimages") 
	{
		int numImg = imgDat.getNumberOfImages();
		if (numImg < 2)
		{
			setErrorString("Images data set must contain at least two images per source");
			return false;
		}

		for (int i = 0 ; i < numImg ; i++)
		{
			int numPoints = imgDat.getNumberOfImagePoints(i);

			if (numPoints < 1)
			{
				setErrorString("An empty image is present");
				return false;
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

	if (typeName != "causticgrid") // ignore the rest
		return true;

	if (m_lastImgIdx < 0)
	{
		setErrorString("No previous extended image information that this caustic grid can refer to");
		return false;
	}

	if (m_lastImgDds != imgDat.getDds() || m_lastImgDs != imgDat.getDs())
	{
		setErrorString("Caustic grid distances are not the same as the ones from the last recorded extended images");
		return false;
	}

	if (!imgDat.hasTriangulation())
	{
		setErrorString("Caustic grid doesn't contain a triangulation");
		return false;
	}

	needCalcDeflections = true;
	needCalcDeflDeriv = true;
	needCalcInverseMag = true;

	int numImages = imgDat.getNumberOfImages();
	int sourcePos = m_sourceIndices.size();

	m_lineSegments.resize(sourcePos+1);
	m_lineSegmentFlags.resize(sourcePos+1);
	m_lineSegmentIntersections.resize(sourcePos+1);
	m_critTriangles.resize(sourcePos+1);

	m_lineSegments[sourcePos].resize(numImages);
	m_lineSegmentFlags[sourcePos].resize(numImages);
	m_lineSegmentIntersections[sourcePos].resize(numImages);
	m_critTriangles[sourcePos].resize(numImages);

	for (int i = 0 ; i < numImages ; i++)
	{
		vector<TriangleIndices> t;

		if (!imgDat.getTriangles(i, t))
		{
			setErrorString("Unable to obtain triangulation data from the caustic grid: " + imgDat.getErrorString());
			return false;
		}

		m_critTriangles[sourcePos][i].resize(t.size());

		int j = 0;
		vector<TriangleIndices>::const_iterator it;

		for (it = t.begin(); it != t.end() ; it++, j++)
		{
			TriangleIndices triangle = *it;

			int lineIndex1 = storeTriangleLine(sourcePos, i, triangle.getIndex(0), triangle.getIndex(1));
			int lineIndex2 = storeTriangleLine(sourcePos, i, triangle.getIndex(1), triangle.getIndex(2));
			int lineIndex3 = storeTriangleLine(sourcePos, i, triangle.getIndex(2), triangle.getIndex(0));
		
			m_critTriangles[sourcePos][i][j] = TriangleIndices(lineIndex1, lineIndex2, lineIndex3);
		}

		m_lineSegmentFlags[sourcePos][i].resize(m_lineSegments[sourcePos][i].size());
		m_lineSegmentIntersections[sourcePos][i].resize(m_lineSegments[sourcePos][i].size());
	}

	m_sourceIndices.push_back(m_lastImgIdx);
	m_critIndices.push_back(idx);

	addImagesDataIndex(m_lastImgIdx);
	addImagesDataIndex(idx);
	addToUsedImagesCount(m_lastImgNumImgs);
	//Note: disabled this again, so we can use the same image for caustics and null space
	//      for the same image
	//m_lastImgIdx = -1; // reset the referenced image
	return true;
}

int FitnessComponent_CausticPenalty::storeTriangleLine(int sourcePos, int imgIndex, int point1, int point2)
{
	int p1 = point1;
	int p2 = point2;

	if (p1 > p2)
	{
		p1 = point2;
		p2 = point1;
	}

	for (int i = 0 ; i < m_lineSegments[sourcePos][imgIndex].size() ; i++)
	{
		if (m_lineSegments[sourcePos][imgIndex][i].first == p1 && m_lineSegments[sourcePos][imgIndex][i].second == p2)
			return i;
	}

	int index = m_lineSegments[sourcePos][imgIndex].size();

	m_lineSegments[sourcePos][imgIndex].push_back(std::pair<int,int>(p1, p2));

	return index;
}

bool FitnessComponent_CausticPenalty::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateCausticPenaltyFitness(iface, m_sourceIndices, m_critIndices,
			                                 m_lineSegments, m_critTriangles, m_lineSegmentFlags,
											 m_lineSegmentIntersections, getCache());
	return true;
}

} // end namespace
