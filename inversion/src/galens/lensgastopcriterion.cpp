
#include "lensgastopcriterion.h"
#include "lensfitnessgeneral.h"
#include "lensgaindividual.h"
#include <sstream>

using namespace std;
using namespace errut;

namespace grale
{

LensGAStopCriterion::LensGAStopCriterion(const std::shared_ptr<LensGAGenomeMutation> &mutation)
	: m_mutation(mutation)
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
	m_fitnessConvergenceFactors = convParams.getConvergenceFactors();
	m_mutationSizes = convParams.getConvergenceSmallMutationSizes();

	if (m_fitnessConvergenceFactors.size() != m_mutationSizes.size() || m_fitnessConvergenceFactors.size() < 1)
		return "Unexpected: invalid convergence or mutation settings (should have been checked before)";
	
	m_pFitnessHistory = make_unique<MultiFitnessHistory>(m_numObjectives, histSize, 0); // Convergence factor will be reset in the next subroutine
	m_convergenceFactorPos = 0;
	updateMutationAndConvergenceInfo();

	return true;
}

void LensGAStopCriterion::updateMutationAndConvergenceInfo()
{
	float mutationSize = m_mutationSizes[m_convergenceFactorPos];
	bool useAbsoluteMutation = (mutationSize < 0)?true:false;

	m_mutation->setAbsoluteMutation(useAbsoluteMutation);
	m_mutation->setMutationAmplitude(mutationSize);
	
	m_pFitnessHistory->reset(m_fitnessConvergenceFactors[m_convergenceFactorPos]);
	m_convergenceFactorPos++;

	if (!useAbsoluteMutation)
		cout << "DEBUG: Setting small mutation with size " << to_string(mutationSize) << endl;
	else
		cout << "DEBUG: Setting large mutations" << endl;
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
			const LensGAFitness &f = static_cast<const LensGAFitness&>(currentBest[j]->fitnessRef());

			assert(i < (int)f.m_fitnesses.size());
			m_pFitnessHistory->processValue(i, f.m_fitnesses[i]);
		}
	}

	std::stringstream ss;
	// Logging generation-1 for comparison with old code
	ss << "Generation " << (generationNumber-1) << ": " << m_pFitnessHistory->getDebugInfo();
	onReport(ss.str());

	if (m_pFitnessHistory->isFinished())
	{
		if (m_convergenceFactorPos == m_fitnessConvergenceFactors.size())
		{
			shouldStop = true;
			m_stopped = true;
		}
		else
			updateMutationAndConvergenceInfo();
	}

	m_pFitnessHistory->advance();

	if (generationNumber > m_maxGenerations)
		shouldStop = true;

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
