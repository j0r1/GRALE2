#include "fitnesscomponent.h"
#include <limits>

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

FitnessComponent::FitnessComponent(const string &name, const shared_ptr<FitnessComponentCache> &pCache) : errut::ErrorBase(name)
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

} // end namespace
