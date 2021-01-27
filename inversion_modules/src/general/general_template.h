#include "lensfitnessgeneral.h"
#include <grale/lensinversiongenome.h>
#include <grale/vector2d.h>
#include <grale/imagesdata.h>
#include <mogal/gafactorymultiobjective.h>
#include <grale/backprojectmatrix.h>
#include <grale/deflectionmatrix.h>
#include <grale/gravitationallens.h>
#include <grale/multifitnesshistory.h>
#include <vector>
#include <memory>
#include <sstream>

using namespace std;

namespace grale
{

template<class Template_ParentClass>
class LensInversionGAFactory_General : public Template_ParentClass, public mogal::GAFactoryMultiObjective
{
public:
#define MAXMUT 0.1f

	LensInversionGAFactory_General() 
	{
		m_numFitness = 0;
		m_pFitnessHistory = 0;
	}

	~LensInversionGAFactory_General()
	{
	}

	LensFitnessObject *createFitnessObject()
	{
		return new LensFitnessGeneral();
	}

	bool subInit(LensFitnessObject *pFitnessObject0)
	{
		auto pFitnessObject = dynamic_cast<const LensFitnessGeneral*>(pFitnessObject0);
		if (!pFitnessObject)
		{
			setErrorString("The lens fitness object is not of expected type");
			return false;
		}

		m_numFitness = pFitnessObject->getNumberOfFitnessComponents();
		// Tell the multi-objective GA how many fitness components there are
		setNumberOfFitnessComponents(m_numFitness);

		int histSize = pFitnessObject->getConvergenceHistorySize();
		m_fitnessConvergenceFactors = pFitnessObject->getConvergenceFactors();
		m_mutationSizes = pFitnessObject->getConvergenceSmallMutationSizes();

		if (m_fitnessConvergenceFactors.size() != m_mutationSizes.size() || m_fitnessConvergenceFactors.size() < 1)
		{
			setErrorString("Unexpected: invalid convergence or mutation settings (should have been checked before)");
			return false;
		}
		
		m_pFitnessHistory = make_unique<MultiFitnessHistory>(m_numFitness, histSize, 0); // Convergence factor will be reset in the next subroutine
		m_convergenceFactorPos = 0;
		updateMutationAndConvergenceInfo();

		return true;
	}

	void updateMutationAndConvergenceInfo()
	{
		m_mutationSize = m_mutationSizes[m_convergenceFactorPos];
		m_useSmallMutation = (m_mutationSize < 0)?false:true;

		m_pFitnessHistory->reset(m_fitnessConvergenceFactors[m_convergenceFactorPos]);
		m_convergenceFactorPos++;

		if (m_useSmallMutation)
			sendMessage("DEBUG: Setting small mutation with size " + to_string(m_mutationSize));
		else
			sendMessage("DEBUG: Setting large mutations");
	}
	
	void onGeneticAlgorithmStep(int generation,  
			                    bool *generationInfoChanged, bool *stopAlgorithm)
	{
		std::vector<mogal::Genome *> bestGenomes;
		getCurrentAlgorithm()->getBestGenomes(bestGenomes);

		for (int i = 0 ; i < m_numFitness ; i++)
		{
			for (int j = 0 ; j < bestGenomes.size() ; j++)
			{
				LensInversionGenome *pGenome = (LensInversionGenome *)bestGenomes[j];

				assert(m_pFitnessHistory);
				m_pFitnessHistory->processValue(i, pGenome->getFitnessValues()[i]);
			}
		}

		assert(m_pFitnessHistory.get());

		std::stringstream ss;
		ss << "Generation " << generation << ": " << m_pFitnessHistory->getDebugInfo();
		sendMessage(ss.str());

		if (m_pFitnessHistory->isFinished())
		{
			if (m_convergenceFactorPos == m_fitnessConvergenceFactors.size())
				*stopAlgorithm = true;
			else
				updateMutationAndConvergenceInfo();
		}

		m_pFitnessHistory->advance();

		if (generation >= Template_ParentClass::getMaximumNumberOfGenerations())
			*stopAlgorithm = true;
	}

	float getChanceMultiplier()
	{
		return 1.0f;
	}
	
	bool useAbsoluteMutation()
	{
		return (m_useSmallMutation)?false:true;
	}
	
	float getMutationAmplitude()
	{
		return m_mutationSize;
	}

	mogal::Genome *selectPreferredGenome(const std::list<mogal::Genome *> &nonDominatedSet) const
	{
		// Respect the ordering of the fitness components
		std::list<mogal::Genome *> genomes = nonDominatedSet;
		std::list<mogal::Genome *> genomes2;
		LensInversionGenome *bestgenome = 0;

		for (int comp = 0 ; comp < m_numFitness ; comp++)
		{
			bestgenome = 0;

			for (auto it = genomes.begin() ; it != genomes.end() ; ++it)
			{
				LensInversionGenome *g = static_cast<LensInversionGenome *>(*it);

				g->setActiveFitnessComponent(comp);
				if (bestgenome == 0)
					bestgenome = g;
				else
				{
					if (g->isFitterThan(bestgenome))
						bestgenome = g;
				}
			}

			genomes2.clear();

			// Ok, now we know a genome that has the lowest 'comp' fitness, let's see it there
			// are others which perfom equally well
			for (auto it = genomes.begin() ; it != genomes.end() ; ++it)
			{
				LensInversionGenome *g = static_cast<LensInversionGenome *>(*it);

				g->setActiveFitnessComponent(comp);
				// if 'bestgenome' is not fitter than 'g' it must have the same fitness with respect to this component
				if (!bestgenome->isFitterThan(g)) 
					genomes2.push_back(g);
			}

			genomes = genomes2;
			genomes2.clear();
		}

		return bestgenome;
	}
protected:
	int m_numFitness;
	float m_mutationSize;
	bool m_useSmallMutation;

	unique_ptr<MultiFitnessHistory> m_pFitnessHistory;
	std::vector<double> m_fitnessConvergenceFactors;
	std::vector<double> m_mutationSizes;
	int m_convergenceFactorPos;
};

} // end namespace
