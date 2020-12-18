#pragma once

#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_PointImagesOverlap : public FitnessComponent
{
public:
	FitnessComponent_PointImagesOverlap(FitnessComponentCache *pCache);
	~FitnessComponent_PointImagesOverlap();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_PointImagesOverlap(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool finalize(double zd, const Cosmology *pCosm) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setScaleType(PointImageScaleType t) { m_scaleType = t; }
private:
	ScaleFactorWorkspace m_workspace;
	std::vector<float> m_scaleFactors;
	PointImageScaleType m_scaleType;
	std::vector<int> m_setScaleGroups, m_useScaleGroups;

	std::map<std::string, std::vector<int>> m_groupnameIndices;
	std::map<int, std::string> m_useScaleNames;
	std::vector<int> m_nogroupnameIndices;

	int m_numSources;
	bool m_isFinalized;
};

class FitnessComponent_PointGroupOverlap : public FitnessComponent
{
public:
	FitnessComponent_PointGroupOverlap(FitnessComponentCache *pCache);
	~FitnessComponent_PointGroupOverlap();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_PointGroupOverlap(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setFitnessRMSType(PointGroupRMSType t) { m_rmsType = t; }
private:
	PointGroupStorage m_pointGroups;
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
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool finalize(double zd, const Cosmology *pCosm) override;
private:
	PointGroupStorage m_pointGroups;
	std::vector<bool> m_rectFlags;
	std::vector<bool> m_useGroupsFlags;
	int m_extendedSourceCount;
};

} // end namespace


