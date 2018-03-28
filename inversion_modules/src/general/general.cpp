/*

  This file is a part of the GRALE Lens inversion modules, which determine
  the fitness of a lens model for the genetic algorithm used by the
  GRALE library.

  Copyright (C) 2008-2012 Jori Liesenborgs

  Contact: jori.liesenborgs@gmail.com
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  
*/

#include "lensfitnessgeneral.h"
#include <grale/gridlensinversiongafactorybase.h>
#include <grale/gridlensinversiongenomebase.h>
#include <grale/vector2d.h>
#include <grale/imagesdata.h>
#include <mogal/gafactorymultiobjective.h>
#include <grale/backprojectmatrixnew.h>
#include <grale/deflectionmatrix.h>
#include <grale/gravitationallens.h>
#include <grale/multifitnesshistory.h>
#include <vector>
#include <sstream>

using namespace std;
using namespace mogal;

namespace grale
{

class GridLensInversionGAFactory_General : public GridLensInversionGAFactoryBase, public GAFactoryMultiObjective
{
public:
#define HISTORYSIZE 250
#define MAXMUT 0.1f

	GridLensInversionGAFactory_General() 
	{
		m_numFitness = 0;
		m_pFitnessHistory = 0;
		m_useSmallMutation = false;
		m_fitnessConvergenceFactors.push_back(0.05);
		m_mutationSizes.push_back(MAXMUT);
		m_convergenceFactorPos = 0;
	}

	~GridLensInversionGAFactory_General()
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
			m_useSmallMutation = false;
			m_mutationSize = MAXMUT;
		}

		vector<Genome *> bestGenomes;
		getCurrentAlgorithm()->getBestGenomes(bestGenomes);

		for (int i = 0 ; i < m_numFitness ; i++)
		{
			for (int j = 0 ; j < bestGenomes.size() ; j++)
			{
				GridLensInversionGenomeBase *pGenome = (GridLensInversionGenomeBase *)bestGenomes[j];

				assert(m_pFitnessHistory);
				m_pFitnessHistory->processValue(i, pGenome->getFitnessValues()[i]);
			}
		}

		assert(m_pFitnessHistory);

		stringstream ss;
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

		if (generation >= getMaximumNumberOfGenerations())
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

	Genome *selectPreferredGenome(const list<Genome *> &nonDominatedSet) const
	{
		// Respect the ordering of the fitness components
		list<Genome *> genomes = nonDominatedSet;
		list<Genome *> genomes2;
		GridLensInversionGenomeBase *bestgenome = 0;

		for (int comp = 0 ; comp < m_numFitness ; comp++)
		{
			bestgenome = 0;

			for (auto it = genomes.begin() ; it != genomes.end() ; ++it)
			{
				GridLensInversionGenomeBase *g = static_cast<GridLensInversionGenomeBase *>(*it);

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
				GridLensInversionGenomeBase *g = static_cast<GridLensInversionGenomeBase *>(*it);

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
	vector<float> m_fitnessConvergenceFactors;
	vector<float> m_mutationSizes;
	int m_convergenceFactorPos;
};

} // end namespace

extern "C"
{
	GAMODS_EXPORT GAFactory *CreateFactoryInstance()
	{
		return new grale::GridLensInversionGAFactory_General();
	}

	GAMODS_EXPORT grale::LensFitnessObject *CreateFitnessObject()
	{
		return new grale::LensFitnessGeneral();
	}
}

