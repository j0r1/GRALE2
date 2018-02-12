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

#include "lensfitnessweak.h"
#include <grale/gridlensinversiongafactorybase.h>
#include <grale/gridlensinversiongenomebase.h>
#include <grale/vector2d.h>
#include <grale/backprojectmatrixnew.h>
#include <grale/deflectionmatrix.h>
#include <grale/gravitationallens.h>
#include <mogal/gafactorysingleobjective.h>
#include <vector>
#include <iostream>

#define DEBUGOUTPUT

namespace grale
{

class GridLensInversionGAFactory_Weak : public GridLensInversionGAFactoryBase, public mogal::GAFactorySingleObjective
{
public:
	GridLensInversionGAFactory_Weak()
	{
		fitnessobject = 0;
		generationcount = 0;
		m_mutationAmplitude = 0.2;
		m_smallestSet = false;
	}

	~GridLensInversionGAFactory_Weak()
	{
		if (fitnessobject)
			delete fitnessobject;
	}

	LensFitnessObject *createFitnessObject()
	{
		return new LensFitnessWeak();
	}

	bool subInit(LensFitnessObject *pFitnessObject)
	{
		return true;
	}

#define HISTORYSIZE 250
	
	void onGeneticAlgorithmStep(int generation, bool *generationInfoChanged, bool *stopAlgorithm)
	{
		if (generation == 0)
		{
			usesmallmutation = false;
			generationcount = 0;
		}

		std::vector<mogal::Genome *> bestGenomes;
		getCurrentAlgorithm()->getBestGenomes(bestGenomes);
		mogal::Genome *overallBest = bestGenomes[0];
		GridLensInversionGenomeBase *g = (GridLensInversionGenomeBase *)overallBest;

		int pos = generation%HISTORYSIZE;

		scorehist[pos] = g->getFitnessValues()[0];

		if (generationcount >= HISTORYSIZE)
		{
			int nextpos = (generation+1)%HISTORYSIZE; // the score on the next position
			                                          // is actually the score of
								  // HISTORYSIZE generations ago

			if (!usesmallmutation)
			{
				if ((scorehist[nextpos]-scorehist[pos])/scorehist[pos] < 1.0f/5.0f)
				{
					usesmallmutation = true;
					m_smallestSet = false;
					m_mutationAmplitude = 0.2;
					generationcount = 0;
				}
			}
			else
			{
				if (!m_smallestSet)
				{
					if ((scorehist[nextpos]-scorehist[pos])/scorehist[pos] < 1.0f/50.0f)
					{
						m_smallestSet = true;
						m_mutationAmplitude = 0.04;
						generationcount = 0;
					}
				}
				else
				{
					if ((scorehist[nextpos]-scorehist[pos])/scorehist[pos] < 1.0f/100.0f)
						*stopAlgorithm = true;
				}
			}
		}

		generationcount++;

		if (generation >= getMaximumNumberOfGenerations())
			*stopAlgorithm = true;

#ifdef DEBUGOUTPUT
		printf("Generation %d: %g\n",generation,(double)scorehist[pos]);
#endif // DEBUGOUTPUT
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
		return m_mutationAmplitude;
	}

protected:
	bool usesmallmutation;
	float scorehist[HISTORYSIZE];
	LensFitnessObject *fitnessobject;
	int generationcount;
	float m_mutationAmplitude;
	bool m_smallestSet;

};

} // end namespace

extern "C"
{
	GAMODS_EXPORT mogal::GAFactory *CreateFactoryInstance()
	{
		return new grale::GridLensInversionGAFactory_Weak();
	}

	GAMODS_EXPORT grale::LensFitnessObject *CreateFitnessObject()
	{
		return new grale::LensFitnessWeak();
	}
}

