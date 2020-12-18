#pragma once

#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_KappaThreshold : public FitnessComponent
{
public:
	FitnessComponent_KappaThreshold(FitnessComponentCache *pCache);
	~FitnessComponent_KappaThreshold();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_KappaThreshold(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	std::vector<float> m_thresholds;
};

class FitnessComponent_KappaGradient : public FitnessComponent
{
public:
	FitnessComponent_KappaGradient(FitnessComponentCache *pCache);
	~FitnessComponent_KappaGradient();
	FitnessComponent *createShortCopy() const override { return new FitnessComponent_KappaGradient(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
private:
	NumericGradientCalculator m_gradCalc;
	int m_srcIdx;
};

} // end namespace

