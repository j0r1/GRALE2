#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_DeflectionAngle : public FitnessComponent
{
public:
	FitnessComponent_DeflectionAngle(FitnessComponentCache *pCache);
	~FitnessComponent_DeflectionAngle();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_DeflectionAngle(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	std::vector<int> m_srcIdx;
};

} // end namespace

