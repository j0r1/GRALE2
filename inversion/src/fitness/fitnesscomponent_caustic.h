#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_CausticPenalty : public FitnessComponent
{
public:
	FitnessComponent_CausticPenalty(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_CausticPenalty();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_CausticPenalty>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &needCalcDeflSecondDeriv) override;
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

