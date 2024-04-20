#pragma once

#include "lensgastopcriterion.h"
#include "inputoutput.h"
#include <eatk/multipopulationevolver.h>
#include <memory>

class Stop : public grale::LensGAStopCriterion
{
public:
	Stop(int popId, size_t generationOffsetForReporting)
		: grale::LensGAStopCriterion(generationOffsetForReporting), m_popId(popId) { }
protected:
	void onReport(const std::string &s)	const override
	{
		if (m_popId < 0)
			WriteLineStdout("GAMESSAGESTR:" + s);
		else
			WriteLineStdout("GAMESSAGESTR: P(" + std::to_string(m_popId) + "):" + s);
	}
private:
	int m_popId;
};

class MultiStop : public eatk::StopCriterion
{
public:
	MultiStop(size_t numEvolvers, size_t generationOffsetForReporting)
	{
		for (int i = 0 ; i < (int)numEvolvers ; i++)
			m_stops.push_back(std::make_shared<Stop>(i, generationOffsetForReporting));
	}

	errut::bool_t initialize(size_t numObjectives, const grale::LensGAConvergenceParameters &convParams)
	{
		errut::bool_t r;
		for (auto &stop : m_stops)
		{
			if (!(r = stop->initialize(numObjectives, convParams)))
				return "Unable to initialize stop criterion for subpopulation: " + r.getErrorString();
		}
		return true;
	}

	errut::bool_t analyze(const eatk::PopulationEvolver &ev, size_t generation, bool &shouldStop)
	{
		if (generation <= 1)
		{
			if (dynamic_cast<const eatk::MultiPopulationEvolver *>(&ev) == nullptr)
				return "Evolver doesn't appear to be a 'MultiPopulationEvolver'";
		}

		const eatk::MultiPopulationEvolver &multiEvolver = static_cast<const eatk::MultiPopulationEvolver &>(ev);
		auto &singleEvolvers = multiEvolver.getSinglePopulationEvolvers();

		if (singleEvolvers.size() != m_stops.size())
			return "Number of single population evolvers (" + std::to_string(singleEvolvers.size()) + ") doesn't match number individual stop criteria (" + std::to_string(m_stops.size()) + ")";

		bool stop = true;
		errut::bool_t r;

		for (size_t i = 0 ; i < m_stops.size() ; i++)
		{
			bool shouldStopSingle = false;

			if (!(r = m_stops[i]->analyze(*singleEvolvers[i], generation, shouldStopSingle)))
				return "Error running stop criterion for population " + std::to_string(i) + ": " + r.getErrorString();
			if (!shouldStopSingle)
				stop = false;
		}

		// Stop if all populations indicate stop
		shouldStop = stop;

		return true;
	}
private:
	std::vector<std::shared_ptr<Stop>> m_stops;
};


