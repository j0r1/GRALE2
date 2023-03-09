#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"
#include "pointgroupstorage.h"

namespace grale
{

class FitnessComponent_ParityPenalty : public FitnessComponent
{
public:
	FitnessComponent_ParityPenalty(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_ParityPenalty();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_ParityPenalty>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	PointGroupStorage m_pointGroups;
};

} // end namespace

