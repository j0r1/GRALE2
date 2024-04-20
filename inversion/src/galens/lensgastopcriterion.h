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
	LensGAStopCriterion(size_t generationNumberOffset = 0);

	errut::bool_t initialize(size_t numObjectives, const LensGAConvergenceParameters &convParams);
	errut::bool_t analyze(const eatk::PopulationEvolver &evolver, size_t generationNumber, bool &shouldStop) override;
protected:
	virtual void onReport(const std::string &s)	const { }
private:
	size_t m_maxGenerations;
	size_t m_numObjectives;
	double m_fitnessConvergenceFactor;
	std::unique_ptr<MultiFitnessHistory> m_pFitnessHistory;
	bool m_stopped;

	std::chrono::time_point<std::chrono::steady_clock> m_lastFitnessReportTime;
	size_t m_generationNumberOffset;
};

}
