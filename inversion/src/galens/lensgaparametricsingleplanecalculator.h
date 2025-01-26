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

class BetaSizeStats;

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
	bool isRandomizingInputPositions() const { return m_randomizingInputPositions; }

	std::shared_ptr<eatk::FitnessComparison> getFitnessComparison() const override;
	
	const std::vector<float> getInitMin() const { return m_initMin; }
	const std::vector<float> getInitMax() const { return m_initMax; }
	const std::vector<float> getHardMin() const { return m_hardMin; }
	const std::vector<float> getHardMax() const { return m_hardMax; }
private:
	errut::bool_t onNewCalculationStart(size_t iteration, size_t genomesForThisCalculator, size_t genomesForPopulationCalculator) override;
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
	bool m_randomizingInputPositions = false;
	std::vector<Vector2Df> m_thetas;
	std::string m_kernelCode;
	std::string m_kernelName;
	size_t m_numGenomesToCalculate = 0;
	std::vector<size_t> m_changeableParamIdx;
	std::vector<int> m_intParams;
	std::vector<float> m_floatParams;
	std::unique_ptr<GravitationalLens> m_templateLens;
	size_t m_numOriginParams = 0;

	std::shared_ptr<std::vector<bool>> m_tracedSourcesFlags;
	std::shared_ptr<std::vector<std::vector<Vector2Df>>> m_tracedSourcesPoints;
	std::vector<std::pair<size_t, size_t>> m_bpPointInfo;

	std::vector<Vector2Df> m_alphas, m_tracedThetas;
	std::vector<float> m_axx, m_ayy, m_axy, m_potential, m_tracedBetaDiffs;
	std::vector<Vector2Df> m_adjustedThetas, m_fullAdjustedThetas;
	std::vector<float> m_betas, m_scaledAlphas;
	std::vector<float> m_scaledAxx, m_scaledAyy, m_scaledAxy;
	std::vector<float> m_scaledPotentials;
	std::vector<double> m_distFrac;
	bool m_needDerivs = false, m_needPotentials = false;
	float m_potScaleConversion = 0;
	bool m_firstCalculationForNewGeneration = false; // Needed to fetch randomized input positions

	std::vector<float> m_initMin, m_initMax;
	std::vector<float> m_hardMin, m_hardMax;
	std::vector<std::shared_ptr<ParameterPrior>> m_priors;
	int m_fitnessToAddPriorTo = -1;

	class ThetaPointMap
	{
	public:
		void clear()
		{
			m_thetaMap.clear();
			m_pointIndex.clear();
		}

		template <class CB1, class CB2>
		errut::bool_t addPoint(Vector2Df pt, bool forceUnique, CB1 newPointCallback, CB2 existingPointCallback)
		{
			errut::bool_t r;
			auto it = m_thetaMap.find(pt);
			if (forceUnique || it == m_thetaMap.end()) // new point
			{
				size_t ptIdx = m_points.size();
				m_points.push_back(pt);
				m_pointIndex.push_back(ptIdx);

				if (!forceUnique) // Don't make future references to unique points
					m_thetaMap[pt] = ptIdx;

				if (!(r = newPointCallback(pt, forceUnique, ptIdx)))
					return r;
			}
			else
			{
				size_t ptIdx = it->second;
				m_pointIndex.push_back(ptIdx);
				if (!(r = existingPointCallback(pt, ptIdx)))
					return r;
			}
			return true;
		}

		const std::vector<size_t> &getPointMapping() const { return m_pointIndex; }
		const std::vector<Vector2Df> &getPoints() const { return m_points; }
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
		std::vector<Vector2Df> m_points;

		std::map<Vector2Df, size_t, Vec2DfComparator> m_thetaMap;
	};

	ThetaPointMap m_pointMap;

	std::unique_ptr<BetaSizeStats> m_stats;
};

} // end namespace
