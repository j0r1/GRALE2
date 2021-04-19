#pragma once

#include "graleconfig.h"
#include "fitnesscomponent.h"

namespace grale
{

class FitnessComponent_WeakLensing : public FitnessComponent
{
public:
	FitnessComponent_WeakLensing(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_WeakLensing();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_WeakLensing>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;

	void setFitnessType(WeakLensingType t) { m_wlType = t; }
private:
	WeakLensingType m_wlType;
	std::vector<float> m_thresholds;
};

class FitnessComponent_WeakLensing_Bayes : public FitnessComponent
{
public:
	FitnessComponent_WeakLensing_Bayes(const std::shared_ptr<FitnessComponentCache> &pCache);
	~FitnessComponent_WeakLensing_Bayes();
	std::unique_ptr<FitnessComponent> createShortCopy() const override { return std::make_unique<FitnessComponent_WeakLensing_Bayes>(nullptr); }

	bool inspectImagesData(int idx, const ImagesDataExtended &imgDat,
			                       bool &needCalcDeflections, bool &needCalcDeflDeriv, bool &needCalcPotential,
			                       bool &needCalcInverseMag, bool &needCalcShear, bool &needCalcConvergence) override;
	bool calculateFitness(const ProjectedImagesInterface &iface, float &fitness) override;
	bool processFitnessOption(const std::string &optionName, const TypedParameter &value) override;
	bool finalize(double zd, const Cosmology *pCosm) override;
private:
	std::vector<int> m_elliptImgs, m_priorDensImages, m_slImages, m_magImages;
	PointGroupStorage m_pointGroups;

	std::vector<std::vector<float>> m_distanceFractionsForZ;

	bool m_redshiftDistributionNeeded;
	float m_howManySigmaFactor;
	int m_numSigmaSamplePoints;

	float m_slSigmaArcsec;
	float m_zLens;
	float m_maxZ; // Keeps track of maximum redshift in input data (uncertainty included, 5sigma)

	DiscreteFunction<float> m_distFracFunction;

	std::vector<float> m_zDistValues;
	std::vector<float> m_zDistMinMax;
	std::unique_ptr<DiscreteFunction<float>> m_zDistFunction;
	float m_zDistSampleCount;
	std::vector<std::pair<float,float>> m_zDistDistFracAndProb;

	std::unique_ptr<DiscreteFunction<float>> m_baDistFunction;
};

} // end namespace

