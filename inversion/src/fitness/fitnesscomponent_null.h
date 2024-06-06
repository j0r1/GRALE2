#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_NullSpacePointImages : public FitnessComponent
{
public:
	FitnessComponent_NullSpacePointImages(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_NullSpacePointImages();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_NullSpacePointImages>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &needCalcDeflSecondDeriv) override;
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
	FitnessComponent_NullSpaceExtendedImages(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_NullSpaceExtendedImages();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_NullSpaceExtendedImages>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &needCalcDeflSecondDeriv) override;
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

} // end namespace

