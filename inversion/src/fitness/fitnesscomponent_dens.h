#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_KappaThreshold : public FitnessComponent
{
public:
	FitnessComponent_KappaThreshold(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_KappaThreshold();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_KappaThreshold>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &needCalcDeflSecondDeriv) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	std::vector<bool> m_minMax;
};

class FitnessComponent_KappaGradient : public FitnessComponent
{
public:
	FitnessComponent_KappaGradient(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_KappaGradient();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_KappaGradient>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence,
								   bool &needCalcDeflSecondDeriv) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	NumericGradientCalculator m_gradCalc;
	int m_srcIdx;
};

} // end namespace

