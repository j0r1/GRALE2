
#include "lensgastopcriterion.h"
#include "lensfitnessgeneral.h"
#include "lensgaindividual.h"
#include <sstream>

using namespace std;
using namespace errut;

namespace grale
{

LensGAStopCriterion::LensGAStopCriterion(size_t generationNumberOffset, bool useGenerationNumberOffsetInStop)
	: m_generationNumberOffset(generationNumberOffset), m_genOffsetInStop(useGenerationNumberOffsetInStop)
{ 
	m_numObjectives = 0; // means not initialized
	m_lastFitnessReportTime = chrono::steady_clock::now();
	m_stopped = false;
}

bool_t LensGAStopCriterion::initialize(size_t numObjectives, const LensGAConvergenceParameters &convParams)
{
	if (numObjectives < 1)
		return "Number of objectives must be at least one";

	if (m_numObjectives > 0)
		return "Already initialized";

	m_numObjectives = numObjectives;
	m_maxGenerations = convParams.getMaximumNumberOfGenerations();

	int histSize = convParams.getConvergenceHistorySize();
	m_fitnessConvergenceFactor = convParams.getConvergenceFactor();

	m_pFitnessHistory = make_unique<MultiFitnessHistory>(m_numObjectives, histSize, m_fitnessConvergenceFactor); // Convergence factor will be reset in the next subroutine

	return true;
}

bool_t LensGAStopCriterion::analyze(const eatk::PopulationEvolver &evolver, size_t generationNumber, bool &shouldStop)
{
	if (m_numObjectives == 0)
		return "Not initialized";

	if (m_stopped) // It's possible that we continue to evolve when we're actually done, in case we're using multiple populations
	{
		shouldStop = true;
		return true;
	}

	assert(m_pFitnessHistory.get());

	auto &currentBest = evolver.getBestIndividuals();

	for (int i = 0 ; i < m_numObjectives ; i++)
	{
		for (int j = 0 ; j < currentBest.size() ; j++)
		{
			const eatk::Fitness &f = static_cast<const eatk::Fitness&>(currentBest[j]->fitnessRef());

			m_pFitnessHistory->processValue(i, f.getRealValue(i));
		}
	}

	std::stringstream ss;
	// Logging generation-1 for comparison with old code
	ss << "Generation " << (generationNumber-1 + m_generationNumberOffset) << ": " << m_pFitnessHistory->getDebugInfo();
	onReport(ss.str());

	if (m_pFitnessHistory->isFinished())
	{
		shouldStop = true;
		m_stopped = true;
	}

	m_pFitnessHistory->advance();

	if (m_genOffsetInStop)
	{
		if (generationNumber + m_generationNumberOffset > m_maxGenerations)
			shouldStop = true;
	}
	else
	{
		if (generationNumber > m_maxGenerations)
			shouldStop = true;
	}

	auto now = chrono::steady_clock::now();
	// For similarity of the logs, in MOGAL a 10 second interval was used
	if (shouldStop || chrono::duration_cast<chrono::milliseconds>(now - m_lastFitnessReportTime).count() > 10000)
	{
		m_lastFitnessReportTime = now;

		// Report ND set		
		stringstream ss;
		ss << "Current best:";
		for (auto &f : currentBest)
			ss << "( " << f->fitness()->toString() << ")";
		onReport(ss.str());
	}

	return true;
}

}
