#ifndef FITNESSCOMPONENT_H

#define FITNESSCOMPONENT_H

#include "pointgroupstorage.h"
#include "fitnessutil.h"
#include <grale/triangleindices.h>
#include <grale/imagesdataextended.h>
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
	virtual FitnessComponent *createShortCopy() const = 0;

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
	virtual bool processFitnessOption(const std::string &optionName, const TypedParameter &value);
	
	const std::vector<int> &getUsedImagesDataIndices() const							{ return m_usedImagesDataIndices; }

	// This will be a tie breaker in case different fitness measures have the same priority
	int getNumberOfUsedImages() const													{ return m_usedImagesCount; }

	std::vector<std::string> getRecognizedTypeNames() const								{ return m_recognizedTypeNames; }
protected:
	void addImagesDataIndex(int idx);
	void addRecognizedTypeName(const std::string &n);
	void addToUsedImagesCount(int c);
	FitnessComponentCache *getCache()													{ return m_pCache; }

	static bool isPointImage(const ImagesData &imgDat);
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
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_PointImagesOverlap(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool finalize() override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setScaleType(PointImageScaleType t) { m_scaleType = t; }
private:
	std::vector<float> m_distanceFractions;
	ScaleFactorWorkspace m_workspace;
	std::vector<float> m_scaleFactors;
	PointImageScaleType m_scaleType;
	std::vector<int> m_setScaleGroups, m_useScaleGroups;

	std::map<std::string, std::vector<int>> m_groupnameIndices;
	std::map<int, std::string> m_useScaleNames;
	std::vector<int> m_nogroupnameIndices;
};

class FitnessComponent_PointGroupOverlap : public FitnessComponent
{
public:
	FitnessComponent_PointGroupOverlap(FitnessComponentCache *pCache);
	~FitnessComponent_PointGroupOverlap();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_PointGroupOverlap(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setFitnessRMSType(PointGroupRMSType t) { m_rmsType = t; }
private:
	PointGroupStorage m_pointGroups;
	std::vector<float> m_distanceFractions;

	PointGroupRMSType m_rmsType;
};

class FitnessComponent_ExtendedImagesOverlap : public FitnessComponent
{
public:
	FitnessComponent_ExtendedImagesOverlap(FitnessComponentCache *pCache);
	~FitnessComponent_ExtendedImagesOverlap();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_ExtendedImagesOverlap(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool finalize() override;
private:
	PointGroupStorage m_pointGroups;
	std::vector<bool> m_rectFlags;
	std::vector<bool> m_useGroupsFlags;
	int m_extendedSourceCount;
};

class FitnessComponent_NullSpacePointImages : public FitnessComponent
{
public:
	FitnessComponent_NullSpacePointImages(FitnessComponentCache *pCache);
	~FitnessComponent_NullSpacePointImages();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_NullSpacePointImages(nullptr); }

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
	ImagesDataExtended m_lastPointGroupImg;

	PointGroupStorage m_pointGroups;
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
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_NullSpaceExtendedImages(nullptr); }

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
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_WeakLensing(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setFitnessType(WeakLensingType t) { m_wlType = t; }
private:
	WeakLensingType m_wlType;
	std::vector<float> m_thresholds;
};

class FitnessComponent_TimeDelay : public FitnessComponent
{
public:
	enum TDFitnessType { Paper2009, NoSrc };

	FitnessComponent_TimeDelay(FitnessComponentCache *pCache);
	~FitnessComponent_TimeDelay();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_TimeDelay(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &storeOrigIntens, bool &storeOrigTimeDelay, bool &storeOrigShear) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setFitnessType(TDFitnessType t) { m_fitnessType = t; }
private:
	TDFitnessType m_fitnessType;
	std::vector<std::pair<int,int>> m_referencePoints; // for relative TD fitness
	std::vector<float> m_tdScaleFactors;
};

class FitnessComponent_KappaThreshold : public FitnessComponent
{
public:
	FitnessComponent_KappaThreshold(FitnessComponentCache *pCache);
	~FitnessComponent_KappaThreshold();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_KappaThreshold(nullptr); }

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
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_CausticPenalty(nullptr); }

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

