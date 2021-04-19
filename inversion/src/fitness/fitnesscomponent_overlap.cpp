#include "fitnesscomponent_overlap.h"
#include "fitnessutil.h"
#include "imagesdataextended.h"
#include "configurationparameters.h"

#include <iostream>

using namespace std;
using namespace errut;

namespace grale
{

// FitnessComponent_PointImagesOverlap

FitnessComponent_PointImagesOverlap::FitnessComponent_PointImagesOverlap(const std::shared_ptr<FitnessComponentCache> &pCache) 
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

FitnessComponent_PointGroupOverlap::FitnessComponent_PointGroupOverlap(const std::shared_ptr<FitnessComponentCache> &pCache) 
	: FitnessComponent("pointgroupoverlap", pCache)
{
	addRecognizedTypeName("pointgroupimages");
	m_rmsType = AllBetas;
}

FitnessComponent_PointGroupOverlap::~FitnessComponent_PointGroupOverlap()
{
}

bool_t FitnessComponent_PointGroupOverlap::extendedOrPointImageDataToPointGroups(const ImagesDataExtended &imgDat,
		                                     PointGroupStorage &pointGroups)
{
	int numImg = imgDat.getNumberOfImages();
	if (numImg < 1)
		return "Images data set doesn't contain any images";
	
	if (imgDat.getNumberOfGroups() > 0)
		pointGroups.add(&imgDat);
	else // no point groups, is this a point image?
	{
		if (!isPointImage(imgDat))
			return "No point groups present, and not a point image";

		string errStr;
		auto grpImg = addGroupsToPointImages(imgDat, errStr);
		if (!grpImg.get())
			return "Couldn't add group info to point image: " + errStr;
		pointGroups.add(grpImg.get());
	}
	return true;
}

bool FitnessComponent_PointGroupOverlap::inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence)
{
	string typeName;

	imgDat.getExtraParameter("type", typeName);
	if (typeName != "pointgroupimages") // ignore
		return true;
	
	bool_t r = extendedOrPointImageDataToPointGroups(imgDat, m_pointGroups);
	if (!r)
	{
		setErrorString(r.getErrorString());
		return false;
	}
	
	// Ok everything checks out

	needCalcDeflections = true;
	needCalcDeflDeriv = true;

	addImagesDataIndex(idx);
	addToUsedImagesCount(imgDat.getNumberOfImages());

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

FitnessComponent_ExtendedImagesOverlap::FitnessComponent_ExtendedImagesOverlap(const std::shared_ptr<FitnessComponentCache> &pCache) 
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

} // end namespace
