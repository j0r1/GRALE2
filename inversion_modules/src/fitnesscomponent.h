#ifndef FITNESSCOMPONENT_H

#define FITNESSCOMPONENT_H

#include "pointgroupstorage.h"
#include <grale/triangleindices.h>
#include <errut/errorbase.h>
#include <grale/polygon2d.h>
#include <vector>

namespace grale
{

class ImagesDataExtended;
class ProjectedImagesInterface;

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
	FitnessComponent(const std::string &name, FitnessComponentCache *pCache);
public:
	~FitnessComponent();

	void setPriority(int p)																{ m_priority = p; }
	int getPriority() const																{ return m_priority; }

	void setPriorityOrder(int idx)														{ m_priorityIndex = idx; }
	int getPriorityOrder() const														{ return m_priorityIndex; }

	virtual bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear);

	// Note that this may be called more than once
	// TODO: prevent this?
	// TODO: check that this is actually the case?
	virtual bool finalize()																{ return true; }
	virtual bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness);
	
	const std::vector<int> &getUsedImagesDataIndices() const							{ return m_usedImagesDataIndices; }

	// This will be a tie breaker in case different fitness measures have the same priority
	int getNumberOfUsedImages() const													{ return m_usedImagesCount; }

	std::vector<std::string> getRecognizedTypeNames() const								{ return m_recognizedTypeNames; }
protected:
	void addImagesDataIndex(int idx);
	void addRecognizedTypeName(const std::string &n);
	void addToUsedImagesCount(int c);
	FitnessComponentCache *getCache()													{ return m_pCache; }
private:
	int m_priority;
	int m_usedImagesCount;
	std::vector<int> m_usedImagesDataIndices;
	std::vector<std::string> m_recognizedTypeNames;
	FitnessComponentCache *m_pCache;
	int m_priorityIndex;
};

class FitnessComponent_PointImagesOverlap : public FitnessComponent
{
public:
	FitnessComponent_PointImagesOverlap(FitnessComponentCache *pCache);
	~FitnessComponent_PointImagesOverlap();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool finalize() override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	std::vector<float> m_distanceFractions;
	std::vector<float> m_workspace;
	std::vector<float> m_scaleFactors;
	std::vector<int> m_setScaleGroups, m_useScaleGroups;

	std::map<std::string, std::vector<int>> m_groupnameIndices;
	std::map<int, std::string> m_useScaleNames;
	std::vector<int> m_nogroupnameIndices;
};

class FitnessComponent_ExtendedImagesOverlap : public FitnessComponent
{
public:
	FitnessComponent_ExtendedImagesOverlap(FitnessComponentCache *pCache);
	~FitnessComponent_ExtendedImagesOverlap();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	PointGroupStorage m_pointGroups;
	std::vector<bool> m_rectFlags;
	std::vector<bool> m_useGroupsFlags;
};

class FitnessComponent_NullSpacePointImages : public FitnessComponent
{
public:
	FitnessComponent_NullSpacePointImages(FitnessComponentCache *pCache);
	~FitnessComponent_NullSpacePointImages();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	int m_lastImgIdx;
	double m_lastImgDds;
	double m_lastImgDs;
	int m_lastImgNumImgs;

	std::vector<int> m_sourceIndices;
	std::vector<int> m_nullIndices;
	std::vector<float> m_nullWeights;
	std::vector<std::vector<TriangleIndices> > m_nullTriangles;
};

// Can still use point images as a penalty where only one image should be
// present
class FitnessComponent_NullSpaceExtendedImages : public FitnessComponent
{
public:
	FitnessComponent_NullSpaceExtendedImages(FitnessComponentCache *pCache);
	~FitnessComponent_NullSpaceExtendedImages();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	int m_lastImgIdx;
	double m_lastImgDds;
	double m_lastImgDs;
	int m_lastImgNumImgs;

	std::vector<int> m_sourceIndices;
	std::vector<int> m_nullIndices;
	std::vector<float> m_nullWeights;
	std::vector<std::vector<TriangleIndices> > m_nullTriangles;
	std::vector<std::vector<double> > m_nullTriangleAreas;
};

class FitnessComponent_WeakLensing : public FitnessComponent
{
public:
	FitnessComponent_WeakLensing(FitnessComponentCache *pCache);
	~FitnessComponent_WeakLensing();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	std::vector<bool> m_reducedShear;
	std::vector<double> m_thresholds;
};

class FitnessComponent_TimeDelay : public FitnessComponent
{
public:
	enum TDFitnessType { Paper2009, ExpI, ExpII };

	FitnessComponent_TimeDelay(FitnessComponentCache *pCache);
	~FitnessComponent_TimeDelay();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;

	void setFitnessType(TDFitnessType t) { m_fitnessType = t; }
	void setUseRelativeDelays(bool f = true) { m_relative = f; }
private:
	TDFitnessType m_fitnessType;
	bool m_relative;
	std::vector<std::pair<int,int>> m_referencePoints; // for relative TD fitness
};

class FitnessComponent_KappaThreshold : public FitnessComponent
{
public:
	FitnessComponent_KappaThreshold(FitnessComponentCache *pCache);
	~FitnessComponent_KappaThreshold();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	std::vector<float> m_thresholds;
};

class FitnessComponent_CausticPenalty : public FitnessComponent
{
public:
	FitnessComponent_CausticPenalty(FitnessComponentCache *pCache);
	~FitnessComponent_CausticPenalty();

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	int storeTriangleLine(int sourcePos, int imgIndex, int point1, int point2);

	int m_lastImgIdx;
	double m_lastImgDds;
	double m_lastImgDs;
	int m_lastImgNumImgs;

	std::vector<int> m_sourceIndices;
	std::vector<int> m_critIndices;
	std::vector<std::vector<std::vector<std::pair<int,int> > > > m_lineSegments;
	std::vector<std::vector<std::vector<bool> > > m_lineSegmentFlags;
	std::vector<std::vector<std::vector<Vector2D<float> > > > m_lineSegmentIntersections;
	std::vector<std::vector<std::vector<TriangleIndices> > > m_critTriangles;
};

} // end namespace

#endif // FITNESSCOMPONENT_H

