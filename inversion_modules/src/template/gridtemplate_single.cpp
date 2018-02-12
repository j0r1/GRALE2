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

// TODO: include what you need
#include "lensfitnessobject.h"
#include <grale/gridlensinversiongafactorybase.h>
#include <grale/gridlensinversiongenomebase.h>
#include <grale/vector2d.h>
#include <mogal/gafactorysingleobjective.h>
#include <vector>
#include <iostream>

namespace grale
{

// A single objective lens inversion module must inherit:
//  - GridLensInversionGAFactoryBase
//  - mogal::GAFactorySingleObjective
class GridLensInversionGAFactory_Template_SingleObjective : public GridLensInversionGAFactoryBase, public mogal::GAFactorySingleObjective
{
public:
	GridLensInversionGAFactory_Template_SingleObjective()
	{
	}

	~GridLensInversionGAFactory_Template_SingleObjective()
	{
	}

	LensFitnessObject *createFitnessObject()
	{
		// TODO: return your own implementation of LensFitnessObject
		//return new LensFitnessTemplate_SingleObjective();
	}

	bool subInit(LensFitnessObject *pFitnessObject)
	{
		// For a single objective algorithm, nothing to be done here
		// For a multi-objective one, you do need to do something 
		// (take a look at the other template)
		return true;
	}


	// These three functions all have to do with the mutations in the genome
	
	// This determines how many basis functions are changed due to the mutation rule.
	// If it's set to 1, on average one basis function is mutated in each genome.
	float getChanceMultiplier()
	{
		return 1.0f;
	}
	
	// If true, the value of the mutated basis function is changed to a random value
	// between zero and one (note that the basis function weights are always rescaled
	// so that they lie between 0 and 1, so this rule corresponds to a new random
	// basis function weight)
	// If false, the current value of the basis function will be changed somewhat.
	// The amount by which it is changed is controlled by 'getMutationAmplitude'
	bool useAbsoluteMutation()
	{
		return true;
	}
	
	// This is only used if the function above 
	float getMutationAmplitude()
	{
		return 0;
	}


	// This function is executed for each generation and must signal the GA to stop at a certain
	// point
	void onGeneticAlgorithmStep(int generation, bool *generationInfoChanged, bool *stopAlgorithm)
	{
		// This is a very simple stopping criterion, a more advanced one
		// should look at how the fitness value changes
		if (generation >= getMaximumNumberOfGenerations())
			*stopAlgorithm = true;

		// Typically this function will also control the return values of 
		// useAbsoluteMutation and getMutationAmplitude, so that the GA at first
		// explores the problem space by using absolute mutations, and then zooms
		// in on a solution by using relative mutations
	}
};

} // end namespace

extern "C"
{
	// A GA module for gravitational lens inversion must have 
	// the functions CreateFactoryInstance and CreateFitnessObject,
	// exposed in a C-like manner.
	
	// This function must create an instance of the class above, which
	// has mogal::GAFactory as a parent.
	GAMODS_EXPORT mogal::GAFactory *CreateFactoryInstance()
	{
		return new grale::GridLensInversionGAFactory_Template_SingleObjective();
	}

	// This function must create a new instance of the LensFitnessObject derived
	// class used by this genetic algorithm. It is used internally in the GA, and
	// also by the 'imgdata/fitness' command in graleshell.
	GAMODS_EXPORT grale::LensFitnessObject *CreateFitnessObject()
	{
		// TODO: return your own implementation of LensFitnessObject
		//return new grale::LensFitnessTemplate_SingleObjective();
	}
}

