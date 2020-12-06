#include "fitnesscomponent.h"
#include "fitnessutil.h"
#include <grale/imagesdataextended.h>
#include <grale/projectedimagesinterface.h>
#include <grale/configurationparameters.h>
#include <grale/cosmology.h>
#include <limits>
#include <list>
#include <sstream>

#include <iostream>

using namespace std;
using namespace errut;

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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	setErrorString("FitnessComponent::inspectImagesData needs to be implemented in a subclass");
	return false;
}

bool FitnessComponent::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	setErrorString("Unrecognized option '" + optionName + "'");
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
	m_numSources = 0;
	m_isFinalized = false;
}

FitnessComponent_PointImagesOverlap::~FitnessComponent_PointImagesOverlap()
{
}

bool FitnessComponent_PointImagesOverlap::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
							   bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
							   bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
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

	if (imgDat.getDds() < 0 || imgDat.getDs() < 0)
	{
		setErrorString("Source/lens and source/observer distances must be positive");
		return false;
	}

	m_numSources++;

	// If there's only a single point, it's not a constraint for this fitness component
	// (but can still be used as a null space constraint)
	if (numImg == 1)
		return true;

	needCalcDeflections = true;

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

bool FitnessComponent_PointImagesOverlap::finalize(double zd, const Cosmology *pCosm)
{
	if (m_isFinalized) // TODO: can this be called more than once?
		return true;

	//cerr << "DEBUG: finalize in FitnessComponent_PointImagesOverlap, for " << (void*)this << endl;

	// Need to build the group indices from the map
	m_setScaleGroups.resize(m_numSources, -1);

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
	m_isFinalized = true;
	return true;
}

bool FitnessComponent_PointImagesOverlap::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "scaletype")
	{
		string valueStr = value.getStringValue();
		if (valueStr == "MinMax")
		{
			setScaleType(MinMax);
			return true;
		}
		if (valueStr == "MAD")
		{
			setScaleType(MAD);
			return true;
		}
		
		setErrorString("Invalid value for '" + optionName + "': '" + valueStr + "'");
		return false;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_PointImagesOverlap::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	//cerr << "In calculateFitness " << (void*)this << endl;
	//cerr << "m_setScaleGroups: " << m_setScaleGroups.size() << endl;
	//cerr << "m_useScaleGroups: " << m_useScaleGroups.size() << endl;
	assert(m_setScaleGroups.size() == m_useScaleGroups.size());

	getScaleFactors_PointImages(iface, getUsedImagesDataIndices(), m_setScaleGroups, m_scaleFactors, m_workspace, m_scaleType);
	fitness = calculateOverlapFitness_PointImages(iface, getUsedImagesDataIndices(), m_useScaleGroups, m_scaleFactors);
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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
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

bool FitnessComponent_PointGroupOverlap::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "rmstype")
	{
		string valueStr = value.getStringValue();

		if (valueStr == "AllBetas")
		{
			setFitnessRMSType(AllBetas);
			return true;
		}

		if (valueStr == "AverageBeta")
		{
			setFitnessRMSType(AverageBeta);
			return true;
		}

		setErrorString("Invalid value for '" + optionName + "': '" + valueStr + "'");
		return false;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_PointGroupOverlap::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateOverlapFitness_PointGroups(m_pointGroups, iface, getUsedImagesDataIndices(),
			                                      m_rmsType);
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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
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

	if (imgDat.getDds() < 0 || imgDat.getDs() < 0)
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

bool FitnessComponent_ExtendedImagesOverlap::finalize(double zd, const Cosmology *pCosm)
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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
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

	if (!imgDat.hasProperty(ImagesData::ShearComponent1) ||
	    !imgDat.hasProperty(ImagesData::ShearComponent2))
	{
		setErrorString("Not all components are present");
		return false;
	}

	if (!imgDat.hasProperty(ImagesData::ShearWeight))
	{
		setErrorString("Shear weights are not present");
		return false;
	}

	needCalcDeflDeriv = true;
	needCalcShear = true;
	needCalcConvergence = true;

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

bool FitnessComponent_WeakLensing::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "type")
	{
		string valueStr = value.getStringValue();

		if (valueStr == "AveragedEllipticities")
		{
			cerr << "Setting WL fitness type to AveragedEllipticities" << endl;
			setFitnessType(AveragedEllipticities);
			return true;
		}
		
		if (valueStr == "RealShear")
		{
			cerr << "Setting WL fitness type to RealShear" << endl;
			setFitnessType(RealShear);
			return true;
		}

		if (valueStr == "RealReducedShear")
		{
			cerr << "Setting WL fitness type to RealReducedShear" << endl;
			setFitnessType(RealReducedShear);
			return true;
		}

		setErrorString("Invalid value for '" + optionName + "': '" + valueStr + "'");
		return false;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_WeakLensing::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateWeakLensingFitness(iface, getUsedImagesDataIndices(), m_wlType, m_thresholds);
	return true;
}

// FitnessComponent_TimeDelay

FitnessComponent_TimeDelay::FitnessComponent_TimeDelay(FitnessComponentCache *pCache) 
	: FitnessComponent("timedelay", pCache),
	  m_fitnessType(NoSrc)
{
	addRecognizedTypeName("pointimages");
	addRecognizedTypeName("extendedimages");
	addRecognizedTypeName("pointgroupimages");
}

FitnessComponent_TimeDelay::~FitnessComponent_TimeDelay()
{
}

bool FitnessComponent_TimeDelay::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
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

	//cerr << "Added TD fitness for " << idx << endl;
	addImagesDataIndex(idx);
	m_tdScaleFactors.push_back(tdScaleFactor);

	return true;
}

bool FitnessComponent_TimeDelay::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "type")
	{
		string valueStr = value.getStringValue();
		if (valueStr == "Paper2009")
		{
			cerr << "Regular Time Delay Fitness (from 2009 article)" << endl;
			setFitnessType(Paper2009);
			return true;
		}
		
		if (valueStr == "NoSrc")
		{
			cerr << "EXPERIMENTAL TIME DELAY FITNESS II - NoSrc" << endl;
			setFitnessType(NoSrc);
			return true;
		}

		setErrorString("Invalid value for '" + optionName + "': '" + valueStr + "'");
		return false;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_TimeDelay::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	if (m_fitnessType == Paper2009)
		fitness = calculateTimeDelayFitnessPaper2009(iface, getUsedImagesDataIndices(), m_tdScaleFactors);
	else if (m_fitnessType == NoSrc)
		fitness = calculateTimeDelayFitnessNoSrc(iface, getUsedImagesDataIndices(), m_tdScaleFactors);
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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
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

	// TODO: use another check? Dds doesn't make much sense in a multi-plane scenario
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

FitnessComponent_KappaGradient::FitnessComponent_KappaGradient(FitnessComponentCache *pCache)
	: FitnessComponent("kappagradient", pCache),
	  m_srcIdx(-1)
{
	addRecognizedTypeName("kappagradientgrids");
}

FitnessComponent_KappaGradient::~FitnessComponent_KappaGradient()
{
}

bool FitnessComponent_KappaGradient::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
								bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
								bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "kappagradientgrids") // ignore
		return true;

	if (m_srcIdx >= 0)
	{
		setErrorString("A grid for the kappa gradient has already been set, only one can be present");
		return false;
	}

	bool_t r;
	if (!(r = m_gradCalc.check(imgDat)))
	{
		setErrorString(r.getErrorString());
		return false;
	}
	
	needCalcDeflDeriv = true; // needed for convergence calculation
	needCalcConvergence = true;

	m_srcIdx = idx;
	addImagesDataIndex(idx);

	return true;
}

bool FitnessComponent_KappaGradient::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	const vector<float> &gradientSizes = m_gradCalc.getGradientSizes(m_srcIdx, iface);
	float sum = 0;

	for (auto x : gradientSizes)
		sum += x;
	
	fitness = sum/gradientSizes.size();
	return true;
}

// FitnessComponent_WeakLensing_Bayes

FitnessComponent_WeakLensing_Bayes::FitnessComponent_WeakLensing_Bayes(FitnessComponentCache *pCache) 
	: FitnessComponent("bayesweaklensing", pCache)
{
	addRecognizedTypeName("bayesellipticities");
	m_redshiftDistributionNeeded = false;
	m_howManySigmaFactor = 3.0f;
	m_numSigmaSamplePoints = 7;
	m_maxZ = 0;
}

FitnessComponent_WeakLensing_Bayes::~FitnessComponent_WeakLensing_Bayes()
{
}

bool FitnessComponent_WeakLensing_Bayes::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "bayesellipticities")
		return true; // ignore

	if (imgDat.getNumberOfImages() != 1)
	{
		setErrorString("Each images data instance can only contain one 'image' for shear info");
		return false;
	}

	if (!imgDat.hasProperty(ImagesData::ShearComponent1) || !imgDat.hasProperty(ImagesData::ShearComponent2))
	{
		setErrorString("No ellipticity info is present");
		return false;
	}

	if (!imgDat.hasProperty(ImagesData::Redshift) || !imgDat.hasProperty(ImagesData::RedshiftUncertainty))
	{
		setErrorString("The points must have redshift and uncertainty settings (z=0 z_sigma=0 is completely unknown, z=X z_sigma=0 is accurate redshift, z=X z_sigma=Y is redshift with uncertainty)");
		return false;
	}

	// Check Dds ==1 and Ds == 1, so that Dds/Ds == 1
	if (imgDat.getDs() != 1 || imgDat.getDds() != 1)
	{
		setErrorString("For this type of image, both Dds and Ds need to be set to exactly 1");
		return false;
	}

	// Check distance fraction
	const int numPoints = imgDat.getNumberOfImagePoints(0);
	if (numPoints == 0)
	{
		setErrorString("No points present");
		return false;
	}

	m_distanceFractionsForZ.push_back(vector<float>());
	vector<float> &distFracForZ = m_distanceFractionsForZ[m_distanceFractionsForZ.size()-1];

	for (int i = 0 ; i < numPoints ; i++)
	{
		// We don't have the lens redshift here, not the cosmology
		// For now just save the redshift if it needs to be converted
		float z = (float)imgDat.getImagePointProperty(ImagesData::Redshift, 0, i);
		float dz = (float)imgDat.getImagePointProperty(ImagesData::RedshiftUncertainty, 0, i);
		if (z < 0 || dz < 0)
		{
			setErrorString("All redshifts and uncertainties must be non-negative (z=" + to_string(z) + " z_sigma=" + to_string(dz) + ")");
			return false;
		}
		
		float dfToCalculate = -1; // negative signals that it does not need to be calculated
		if (z == 0) // indicates unknown distance fraction, will need probability info
		{
			if (dz != 0)
			{
				setErrorString("For unknown redshifts, the uncertainty must be set to zero (z_sigma=" + to_string(dz) + ")");
				return false;
			}
			m_redshiftDistributionNeeded = true;
		}
		else
		{
			// Keep track of maximum redshift
			// TODO: is 5 sigma a good upper limit to add?
			m_maxZ = std::max(m_maxZ, z + 5.0f*dz);
			if (dz == 0)
				dfToCalculate = z;
		}

		distFracForZ.push_back(dfToCalculate);
	}
	
	needCalcDeflDeriv = true;
	// These two are not needed, they are calculated from the derivatives
	// in the fitness function
	// needCalcShear = true;
	// needCalcConvergence = true;

	addImagesDataIndex(idx);

	return true;
}

// This is called before inspect
bool FitnessComponent_WeakLensing_Bayes::processFitnessOption(const std::string &optionName, const TypedParameter &value)
{
	if (optionName == "sigmafactor")
	{
		if (value.isArray() || !value.isReal())
		{
			setErrorString("This should be a real number");
			return false;
		}
		m_howManySigmaFactor = value.getRealValue();
		return true;
	}

	if (optionName == "sigmasteps")
	{
		if (value.isArray() || !value.isInteger())
		{
			setErrorString("This should be an integer");
			return false;
		}

		m_numSigmaSamplePoints = value.getIntegerValue();
		if (m_numSigmaSamplePoints < 2)
		{
			setErrorString("At least two points are needed");
			return false;
		}
		return true;
	}

	// TODO: add b/a distribution, check between 0 and 1, calculate area and normalize
	if (optionName == "b_over_a_distribution")
	{
		if (value.isEmpty()) // Allow empty default, check later if it's been set if needed
			return true;

		if (!value.isArray() || !value.isReal() || value.getNumberOfEntries() < 2)
		{
			setErrorString("This should either be an array with at least two real numbers");
			return false;
		}

		// Normalize it ; b/a goes from 0 to 1
		const vector<double> &probdist = value.getRealValues();
		double area = 0;
		for (int i = 0 ; i < probdist.size()-1 ; i++)
		{
			double x0 = (double)i/(double)(probdist.size()-1);
			double x1 = (double)(i+1)/(double)(probdist.size()-1);
			double y0 = probdist[i];
			double y1 = probdist[i+1];
			if (y0 < 0 || y1 < 0)
			{
				setErrorString("All values in b/a distribution must be non-negative");
				return false;
			}

			area += 0.5*(y0+y1)*(x1-x0);
		}
		vector<float> floatProbdist;
		for (auto p : probdist)
			floatProbdist.push_back((float)(p / area));

		m_baDistFunction = make_unique<DiscreteFunction<float>>();
		auto r = m_baDistFunction->init(0.0f, 1.0f, floatProbdist);
		if (!r)
		{
			setErrorString("Unexpected: unable to initialize b/a prob dist function: " + r.getErrorString());
			return false;
		}

		return true;
	}

	// TODO: this should become redshift distribution
	if (optionName == "distfracdistribution")
	{
		if (value.isEmpty()) // Ok, nothing set, but check later if we need it
			return true;
		
		if (!value.isArray() || !value.isReal())
		{
			setErrorString("This should either be empty or an array of real numbers");
			return false;
		}

		const vector<double> &values = value.getRealValues();
		if (values.size() < 1 || values.size()%2 != 0)
		{
			setErrorString("An even number of values should be present, corresponding to pairs of a Dds/Ds fraction and a (positive) weight");
			return false;
		}

		vector<pair<double,double>> weights;
		double weightSum = 0;
		for (size_t i = 0 ; i < values.size() ; i += 2)
		{
			double dfrac = values[i];
			double w = values[i+1];
			if (dfrac < 0 || dfrac > 1.0)
			{
				setErrorString("All distance fractions must be positive (detected " + to_string(dfrac) + ")");
				return false;
			}

			if (w <= 0)
			{
				setErrorString("All weights must be strictly positive (detected " + to_string(w) + ")");
				return false;
			}

			weights.push_back({dfrac, w});
			weightSum += w;
		}

		// Make sure the weights are normalized
		for (auto &dw : weights)
			dw.second /= weightSum;

		m_distanceFractionWeights.clear();
		for (auto dw : weights)
			m_distanceFractionWeights.push_back({(float)dw.first, (float)dw.second});
		return true;
	}

	setErrorString("Unknown option");
	return false;
}

bool FitnessComponent_WeakLensing_Bayes::calculateFitness(const ProjectedImagesInterface &iface, float &fitness)
{
	fitness = calculateWeakLensingFitness_Bayes(iface, getUsedImagesDataIndices(), m_distanceFractionWeights,
	                                            m_howManySigmaFactor, m_numSigmaSamplePoints);
	return true;
}

bool FitnessComponent_WeakLensing_Bayes::finalize(double zd, const Cosmology *pCosm)
{
	if (pCosm == nullptr)
	{
		setErrorString("Cosmological model needs to be set for this component");
		return false;
	}

	if (zd <= 0) // Should not happen of course, just for safety an extra check
	{
		setErrorString("Lens redshift should be positive");
		return false;
	}

	// Determine a discrete distance fraction function
	if (m_maxZ < zd)
		m_maxZ = zd*1.1f; // TODO: should probably not happen anyway
	int distFracPoints = 8192; // TODO: what's a good value here? get this from config option?
	vector<float> distFracs(distFracPoints);
	for (int i = 0 ; i < distFracPoints ; i++)
	{
		float frac = (float)i/(float)(distFracPoints-1);
		float zs = (1.0f-frac)*zd + frac*m_maxZ;
		distFracs[i] = (float)(pCosm->getAngularDiameterDistance(zd, zs)/pCosm->getAngularDiameterDistance(zs));
	}

	auto r = m_distFracFunction.init(zd, m_maxZ, distFracs);
	if (!r)
	{
		setErrorString("Unexpected: couldn't init discrete distance fraction function: " + r.getErrorString());
		return false;
	}

	// Pre-calculate Dds/Ds for points where it's useful
	for (auto &sources : m_distanceFractionsForZ)
	{
		for (auto &zsToDf : sources)
		{
			// If the value is positive, it can be pre-calculated to a distance fraction
			if (zsToDf > 0)
				zsToDf = m_distFracFunction(zsToDf); // TODO: use exact calculation here instead?
		}
	}

	if (m_redshiftDistributionNeeded)
	{
		// TODO: this should become the redshift distribution
		if (m_distanceFractionWeights.size() == 0)
		{
			setErrorString("No distance fraction distribution is set, but not all distance fractions (Dds/Ds) are known");
			return false;
		}
	}

	if (m_baDistFunction.get() == nullptr)
	{
		setErrorString("No b/a distribution has been set");
		return false;
	}


	return true;
}

} // end namespace
