#pragma once

#include "graleconfig.h"
#include "lensgagenomecalculator.h"
#include "lensfitnessobject.h"
#include "lensinversionparametersparametricsingleplane.h"
#include "oclcalculatedbackprojector.h"
#include <cassert>
#include <map>

namespace grale
{

// TODO: for now, all properties (deflection, derivatives, potential) are calculated
//       for all image points
class LensGAParametricSinglePlaneCalculator : public LensGAGenomeCalculator
{
public:
	LensGAParametricSinglePlaneCalculator(std::unique_ptr<LensFitnessObject> fitObj);
	~LensGAParametricSinglePlaneCalculator();

	const LensFitnessObject &getFitnessObject() const { return *m_fitObj; }

	errut::bool_t init(const LensInversionParametersBase &params) override;
	errut::bool_t createLens(const eatk::Genome &genome, std::unique_ptr<GravitationalLens> &lens) const override;
	size_t getNumberOfObjectives() const override { assert(m_fitObj.get()); return m_fitObj->getNumberOfFitnessComponents(); }

	std::shared_ptr<eatk::FitnessComparison> getFitnessComparison() const override;
	
	const std::vector<float> getInitMin() const { return m_initMin; }
	const std::vector<float> getInitMax() const { return m_initMax; }
	const std::vector<float> getHardMin() const { return m_hardMin; }
	const std::vector<float> getHardMax() const { return m_hardMax; }
private:
	errut::bool_t onNewCalculationStart(size_t genomesForThisCalculator, size_t genomesForPopulationCalculator) override;
	errut::bool_t startNewCalculation(const eatk::Genome &genome0) override;
	errut::bool_t pollCalculate(const eatk::Genome &genome0, eatk::Fitness &fitness0) override;

	std::unique_ptr<LensFitnessObject> m_fitObj;
	size_t m_numObjectives = 0;

	std::unique_ptr<OclCalculatedBackProjector> m_oclBp;
	bool m_init = false;
	bool m_infOnBoundsViolation = false;
	double m_angularScale = 0;
	double m_potScale = 0;
	int m_devIdx = 0;
	std::vector<Vector2Df> m_thetas;
	std::string m_kernelCode;
	std::string m_kernelName;
	size_t m_numGenomesToCalculate = 0;
	std::vector<size_t> m_changeableParamIdx;
	std::vector<int> m_intParams;
	std::vector<float> m_floatParams;
	std::unique_ptr<GravitationalLens> m_templateLens;

	std::vector<Vector2Df> m_alphas;
	std::vector<float> m_axx, m_ayy, m_axy, m_potential;
	std::vector<float> m_betas, m_scaledAlphas;
	std::vector<float> m_scaledAxx, m_scaledAyy, m_scaledAxy;
	std::vector<float> m_scaledPotentials;
	std::vector<double> m_distFrac;
	bool m_needDerivs = false, m_needPotentials = false;
	float m_potScaleConversion = 0;

	std::vector<float> m_initMin, m_initMax;
	std::vector<float> m_hardMin, m_hardMax;

	class ThetaPointMap
	{
	public:
		void clear()
		{
			m_thetaMap.clear();
			m_pointIndex.clear();
		}

		bool addPoint(Vector2Df pt)
		{
			auto it = m_thetaMap.find(pt);
			if (it == m_thetaMap.end()) // new point
			{
				size_t ptIdx = m_thetaMap.size();
				m_pointIndex.push_back(ptIdx);

				m_thetaMap[pt] = ptIdx;
				return true;
			}
			else
			{
				size_t ptIdx = it->second;
				m_pointIndex.push_back(ptIdx);
				return false;
			}
		}

		const std::vector<size_t> &getPointMapping() const { return m_pointIndex; }
	private:
		struct Vec2DfComparator
		{
		    bool operator()(const Vector2Df a, const Vector2Df b) const
			{
				if (a.getX() == b.getX())
					return a.getY() < b.getY();
				return a.getX() < b.getX();
		    }
		};
		std::vector<size_t> m_pointIndex;
		std::map<Vector2Df, size_t, Vec2DfComparator> m_thetaMap;
	};

	ThetaPointMap m_pointMap;
};

} // end namespace
