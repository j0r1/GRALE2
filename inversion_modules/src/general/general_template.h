#include "lensfitnessgeneral.h"
#include <grale/lensinversiongenome.h>
#include <grale/vector2d.h>
#include <grale/imagesdata.h>
#include <mogal/gafactorymultiobjective.h>
#include <grale/backprojectmatrixnew.h>
#include <grale/deflectionmatrix.h>
#include <grale/gravitationallens.h>
#include <grale/multifitnesshistory.h>
#include <vector>
#include <sstream>

namespace grale
{

template<class Template_ParentClass>
class LensInversionGAFactory_General : public Template_ParentClass, public mogal::GAFactoryMultiObjective
{
public:
#define HISTORYSIZE 250
#define MAXMUT 0.1f

	LensInversionGAFactory_General() 
	{
		m_numFitness = 0;
		m_pFitnessHistory = 0;
		m_useSmallMutation = false;
		m_fitnessConvergenceFactors.push_back(0.05);
		m_mutationSizes.push_back(MAXMUT);
		m_convergenceFactorPos = 0;
	}

	~LensInversionGAFactory_General()
	{
		delete m_pFitnessHistory;

	}

	LensFitnessObject *createFitnessObject()
	{
		return new LensFitnessGeneral();
	}

	bool subInit(LensFitnessObject *pFitnessObject)
	{
		m_numFitness = pFitnessObject->getNumberOfFitnessComponents();

		delete m_pFitnessHistory;
		m_pFitnessHistory = new MultiFitnessHistory(m_numFitness, HISTORYSIZE, 0.1);

		// Tell the multi-objective GA how many fitness components there are
		setNumberOfFitnessComponents(m_numFitness);
		return true;
	}
	
	void onGeneticAlgorithmStep(int generation,  
			                    bool *generationInfoChanged, bool *stopAlgorithm)
	{
		if (generation == 0)
		{
			if (getenv("GRALE_DEBUG_INSTANTSMALLMUTATION"))
			{
				m_useSmallMutation = true;
				sendMessage("DEBUG: Forcing small mutation from start");
			}
			else
				m_useSmallMutation = false;
	
			m_mutationSize = MAXMUT;
		}

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

		assert(m_pFitnessHistory);

		std::stringstream ss;
		ss << "Generation " << generation << ": " << m_pFitnessHistory->getDebugInfo();
		sendMessage(ss.str());

		if (m_pFitnessHistory->isFinished())
		{
			if (m_convergenceFactorPos == m_fitnessConvergenceFactors.size())
				*stopAlgorithm = true;
			else
			{
				m_useSmallMutation = true;
				m_mutationSize = m_mutationSizes[m_convergenceFactorPos];
				m_pFitnessHistory->reset(m_fitnessConvergenceFactors[m_convergenceFactorPos]);
				m_convergenceFactorPos++;
			}
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

	MultiFitnessHistory *m_pFitnessHistory;
	std::vector<float> m_fitnessConvergenceFactors;
	std::vector<float> m_mutationSizes;
	int m_convergenceFactorPos;
};

} // end namespace
