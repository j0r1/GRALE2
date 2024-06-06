#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_TimeDelay : public FitnessComponent
{
public:
	enum TDFitnessType { Paper2009, NoSrc };

	FitnessComponent_TimeDelay(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_TimeDelay();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_TimeDelay>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &needCalcDeflSecondDeriv) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setFitnessType(TDFitnessType t) { m_fitnessType = t; }
	bool finalize(double zd, const Cosmology *pCosm) override;
private:
	TDFitnessType m_fitnessType;
	std::vector<std::pair<int,int>> m_referencePoints; // for relative TD fitness
	std::vector<float> m_tdScaleFactors;
	float m_nosrcCutoff; // If the NoSrc fitness value is below this, it is set to zero (can help multi-objective)
};

} // end namespace

