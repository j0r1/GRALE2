#pragma once

#include "graleconfig.h"
#include "lensgagenomemutation.h"
#include "multifitnesshistory.h"
#include "lensgaconvergenceparameters.h"
#include <eatk/stopcriterion.h>
#include <chrono>

namespace grale
{

class LensFitnessObject;

class LensGAStopCriterion : public eatk::StopCriterion
{
public:
	LensGAStopCriterion(const std::shared_ptr<LensGAGenomeMutation> &mutation);

	errut::bool_t initialize(size_t numObjectives, const LensGAConvergenceParameters &convParams);
	errut::bool_t analyze(const std::vector<std::shared_ptr<eatk::Individual>> &currentBest, size_t generationNumber, bool &shouldStop) override;
protected:
	virtual void onReport(const std::string &s)	const { }
private:
	void updateMutationAndConvergenceInfo();

	size_t m_maxGenerations;
	std::shared_ptr<LensGAGenomeMutation> m_mutation;
	size_t m_numObjectives;
	std::vector<double> m_fitnessConvergenceFactors;
	std::vector<double> m_mutationSizes;
	std::unique_ptr<MultiFitnessHistory> m_pFitnessHistory;
	int m_convergenceFactorPos;

	std::chrono::time_point<std::chrono::steady_clock> m_lastFitnessReportTime;
};

}
