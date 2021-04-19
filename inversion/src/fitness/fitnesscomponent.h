#pragma once

#include "graleconfig.h"
#include "pointgroupstorage.h"
#include "fitnessutil.h"
#include "numericgradientcalculator.h"
#include "triangleindices.h"
#include "imagesdataextended.h"
#include <errut/errorbase.h>
#include "polygon2d.h"
#include "discretefunction.h"
#include <vector>
#include <memory>

namespace grale
{

class ImagesDataExtended;
class ProjectedImagesInterface;
class Cosmology;

class FitnessComponentCache
{
public:
	FitnessComponentCache(int num);
	~FitnessComponentCache();

	void clear();

	void setEstimatedSource(int idx, const Polygon2D<float> &poly, float minX, float maxX, float minY, float maxY);
	bool getEstimatedSource(int idx, Polygon2D<float> &poly, float &minX, float &maxX, float &minY, float &maxY);

	void setEstimatedSourceScale(int idx, float scale);
	bool getEstimatedSourceScale(int idx, float &scale);
private:
	struct SourceShape
	{
		SourceShape() { m_isEstSrcSet = false; }

		bool m_isEstSrcSet;
		Polygon2D<float> m_estimatedSource;
		float m_srcMinX, m_srcMaxX;
		float m_srcMinY, m_srcMaxY;
	};

	struct SourceScale
	{
		SourceScale() { m_isSet = false; }

		bool m_isSet;
		float m_scale;
	};

	bool m_clear;
	int m_numImgDatas;
	std::vector<SourceShape> m_estimatedSourceShapes;
	std::vector<SourceScale> m_sourceScales;
};

class FitnessComponent : public errut::ErrorBase
{
protected:
	FitnessComponent(const std::string &name, const std::shared_ptr<FitnessComponentCache> &pCache);
public:
	~FitnessComponent();
	virtual std::unique_ptr<FitnessComponent> createShortCopy() const = 0;

	void setPriority(int p)																{ m_priority = p; }
	int getPriority() const																{ return m_priority; }

	void setPriorityOrder(int idx)														{ m_priorityIndex = idx; }
	int getPriorityOrder() const														{ return m_priorityIndex; }

	virtual bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence);

	// Note that this may be called more than once
	// TODO: prevent this?
	// TODO: check that this is actually the case?
	// Note that cosmology may be nullptr
	virtual bool finalize(double zd, const Cosmology *pCosm)							{ return true; }
	virtual bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness);
	virtual bool processFitnessOption(const std::string &optionName, const TypedParameter &value);
	
	const std::vector<int> &getUsedImagesDataIndices() const							{ return m_usedImagesDataIndices; }

	// This will be a tie breaker in case different fitness measures have the same priority
	int getNumberOfUsedImages() const													{ return m_usedImagesCount; }

	std::vector<std::string> getRecognizedTypeNames() const								{ return m_recognizedTypeNames; }

	static bool isPointImage(const ImagesData &imgDat);
protected:
	void addImagesDataIndex(int idx);
	void addRecognizedTypeName(const std::string &n);
	void addToUsedImagesCount(int c);
	FitnessComponentCache *getCache()													{ return m_pCache.get(); }
private:
	int m_priority;
	int m_usedImagesCount;
	std::vector<int> m_usedImagesDataIndices;
	std::vector<std::string> m_recognizedTypeNames;
	std::shared_ptr<FitnessComponentCache> m_pCache;
	int m_priorityIndex;
};

} // end namespace

