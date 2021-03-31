
#include "lensgastopcriterion.h"
#include "lensfitnessgeneral.h"
#include "lensgaindividual.h"
#include <sstream>

using namespace std;
using namespace errut;

namespace grale
{

LensGAStropCriterion::LensGAStropCriterion(size_t maxGenerations,
						const std::shared_ptr<LensGAGenomeMutation> &mutation)
	: m_maxGenerations(maxGenerations), m_mutation(mutation)
{ 
	m_numObjectives = 0; // means not initialized
}

bool_t LensGAStropCriterion::initialize(const LensFitnessObject &fitnessObject)
{
	if (m_numObjectives > 0)
		return "Already initialized";

	auto pFitnessObject = dynamic_cast<const LensFitnessGeneral*>(&fitnessObject);
	if (!pFitnessObject)
		return "Not a LensFitnessGeneral object";

	m_numObjectives = fitnessObject.getNumberOfFitnessComponents();

	int histSize = pFitnessObject->getConvergenceHistorySize();
	m_fitnessConvergenceFactors = pFitnessObject->getConvergenceFactors();
	m_mutationSizes = pFitnessObject->getConvergenceSmallMutationSizes();

	if (m_fitnessConvergenceFactors.size() != m_mutationSizes.size() || m_fitnessConvergenceFactors.size() < 1)
		return "Unexpected: invalid convergence or mutation settings (should have been checked before)";
	
	m_pFitnessHistory = make_unique<MultiFitnessHistory>(m_numObjectives, histSize, 0); // Convergence factor will be reset in the next subroutine
	m_convergenceFactorPos = 0;
	updateMutationAndConvergenceInfo();

	return true;
}

void LensGAStropCriterion::updateMutationAndConvergenceInfo()
{
	float mutationSize = m_mutationSizes[m_convergenceFactorPos];
	bool useAbsoluteMutation = (mutationSize < 0)?true:false;

	m_mutation->setAbsoluteMutation(useAbsoluteMutation);
	m_mutation->setMutationAmplitude(mutationSize);
	
	m_pFitnessHistory->reset(m_fitnessConvergenceFactors[m_convergenceFactorPos]);
	m_convergenceFactorPos++;

	if (!useAbsoluteMutation)
		cerr << "DEBUG: Setting small mutation with size " << to_string(mutationSize) << endl;
	else
		cerr << "DEBUG: Setting large mutations" << endl;
}

bool_t LensGAStropCriterion::analyze(const std::vector<std::shared_ptr<eatk::Individual>> &currentBest, size_t generationNumber, bool &shouldStop)
{
	if (m_numObjectives == 0)
		return "Not initialized";

	assert(m_pFitnessHistory.get());

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
	cerr << ss.str() << endl;

	if (m_pFitnessHistory->isFinished())
	{
		if (m_convergenceFactorPos == m_fitnessConvergenceFactors.size())
			shouldStop = true;
		else
			updateMutationAndConvergenceInfo();
	}

	m_pFitnessHistory->advance();

	if (generationNumber > m_maxGenerations)
		shouldStop = true;

	return true;
}

}