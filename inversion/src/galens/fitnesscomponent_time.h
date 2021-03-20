#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_TimeDelay : public FitnessComponent
{
public:
	enum TDFitnessType { Paper2009, NoSrc };

	FitnessComponent_TimeDelay(FitnessComponentCache *pCache);
	~FitnessComponent_TimeDelay();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_TimeDelay(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setFitnessType(TDFitnessType t) { m_fitnessType = t; }
	bool finalize(double zd, const Cosmology *pCosm) override;
private:
	TDFitnessType m_fitnessType;
	std::vector<std::pair<int,int>> m_referencePoints; // for relative TD fitness
	std::vector<float> m_tdScaleFactors;
};

} // end namespace

