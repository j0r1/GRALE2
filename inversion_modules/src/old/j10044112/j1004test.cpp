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

#include "lensfitnessj1004.h"
#include <grale/lensinversiongafactorysingleplanecpu.h>
#include <grale/lensinversiongenome.h>
#include <grale/vector2d.h>
#include <grale/imagesdata.h>
#include <mogal/gafactorymultiobjective.h>
#include <grale/backprojectmatrix.h>
#include <grale/deflectionmatrix.h>
#include <grale/gravitationallens.h>
#include <grale/multifitnesshistory.h>
#include <vector>
#include <iostream>

namespace grale
{

class LensInversionGAFactorySinglePlaneCPU_J1004Test : public LensInversionGAFactorySinglePlaneCPU, public mogal::GAFactoryMultiObjective
{
public:
#define HISTORYSIZE 150
#define MAXMUT 0.1f

	LensInversionGAFactorySinglePlaneCPU_J1004Test() : m_fitnessHistory(4, HISTORYSIZE, 0.1)
	{
		usesmallmutation = false;
		m_fitnessConvergenceFactors.push_back(0.05);
		m_mutationSizes.push_back(MAXMUT);
		m_convergenceFactorPos = 0;
	}

	~LensInversionGAFactorySinglePlaneCPU_J1004Test()
	{
	}

	LensFitnessObject *createFitnessObject()
	{
		return new LensFitnessJ1004Test();
	}

	bool subInit(LensFitnessObject *pFitnessObject)
	{
		// Tell the multi-objective GA how many fitness components there are
		setNumberOfFitnessComponents(pFitnessObject->getNumberOfFitnessComponents());
		return true;
	}
	
	void onGeneticAlgorithmStep(int generation,  
			                    bool *generationInfoChanged, bool *stopAlgorithm)
	{
		if (generation == 0)
		{
			usesmallmutation = false;
			mutationsize = MAXMUT;
			startgeneration = 0;
		}

		std::vector<mogal::Genome *> bestGenomes;
		getCurrentAlgorithm()->getBestGenomes(bestGenomes);

		for (int i = 0 ; i < 4 ; i++)
		{
			for (int j = 0 ; j < bestGenomes.size() ; j++)
			{
				LensInversionGenome *pGenome = (LensInversionGenome *)bestGenomes[j];

				m_fitnessHistory.processValue(i, pGenome->getFitnessValues()[i]);
			}
		}

		std::cerr << "Generation " << generation << ":";
		m_fitnessHistory.printDebugInfo();

		if (m_fitnessHistory.isFinished())
		{
			if (m_convergenceFactorPos == m_fitnessConvergenceFactors.size())
				*stopAlgorithm = true;
			else
			{
				usesmallmutation = true;
				mutationsize = m_mutationSizes[m_convergenceFactorPos];
				m_fitnessHistory.reset(m_fitnessConvergenceFactors[m_convergenceFactorPos]);
				m_convergenceFactorPos++;
			}
			//*stopAlgorithm = true;
		}

		m_fitnessHistory.advance();

		if (generation >= getMaximumNumberOfGenerations())
			*stopAlgorithm = true;
	}

	float getChanceMultiplier()
	{
		return 1.0f;
	}
	
	bool useAbsoluteMutation()
	{
		return (usesmallmutation)?false:true;
	}
	
	float getMutationAmplitude()
	{
		return mutationsize;
	}

	mogal::Genome *selectPreferredGenome(const std::list<mogal::Genome *> &nonDominatedSet) const
	{
		// First select the genomes with the lowest critical line crossings
		std::list<mogal::Genome *> genomes, genomes2;
		std::list<mogal::Genome *>::const_iterator it;
		grale::LensInversionGenome *bestgenome = 0;

		for (it = nonDominatedSet.begin() ; it != nonDominatedSet.end() ; ++it)
		{
			grale::LensInversionGenome *g = (grale::LensInversionGenome *)(*it);

			g->setActiveFitnessComponent(1);
			if (bestgenome == 0)
				bestgenome = g;
			else
			{
				if (g->isFitterThan(bestgenome))
					bestgenome = g;
			}
		}

		if (bestgenome != 0)
		{
			for (it = nonDominatedSet.begin() ; it != nonDominatedSet.end() ; ++it)
			{
				grale::LensInversionGenome *g = (grale::LensInversionGenome *)(*it);

				if (!bestgenome->isFitterThan(g))
					genomes.push_back(g);
			}
		}

		// Then select the ones with the lowest null space value
		bestgenome = 0;
		for (it = genomes.begin() ; it != genomes.end() ; ++it)
		{
			grale::LensInversionGenome *g = (grale::LensInversionGenome *)(*it);

			g->setActiveFitnessComponent(2);
			if (bestgenome == 0)
				bestgenome = g;
			else
			{
				if (g->isFitterThan(bestgenome))
					bestgenome = g;
			}
		}

		if (bestgenome != 0)
		{
			for (it = genomes.begin() ; it != genomes.end() ; ++it)
			{
				grale::LensInversionGenome *g = (grale::LensInversionGenome *)(*it);

				if (!bestgenome->isFitterThan(g))
					genomes2.push_back(g);
			}
		}

		// Then select the one with the lowest overlap value
		bestgenome = 0;
		for (it = genomes2.begin() ; it != genomes2.end() ; ++it)
		{
			grale::LensInversionGenome *g = (grale::LensInversionGenome *)(*it);

			g->setActiveFitnessComponent(0);
			if (bestgenome == 0)
				bestgenome = g;
			else
			{
				if (g->isFitterThan(bestgenome))
					bestgenome = g;
			}
		}

		return bestgenome;
	}
protected:
	float mutationsize;
	bool usesmallmutation;
	LensFitnessObject *fitnessobject;
	int startgeneration;

	MultiFitnessHistory m_fitnessHistory;
	std::vector<float> m_fitnessConvergenceFactors;
	std::vector<float> m_mutationSizes;
	int m_convergenceFactorPos;
};

} // end namespace

extern "C"
{
	GAMODS_EXPORT mogal::GAFactory *CreateFactoryInstance()
	{
		return new grale::LensInversionGAFactorySinglePlaneCPU_J1004Test();
	}

	GAMODS_EXPORT grale::LensFitnessObject *CreateFitnessObject()
	{
		return new grale::LensFitnessJ1004Test();
	}
}

