/*

  This file is a part of GRALE, a library to facilitate the simulation
  and inversion of gravitational lenses.

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

#include "gridlensinversiongenomebase.h"
#include "lensinversiongafactorycommon.h"
#include "simpleuniformdistribution.h"
#include "constants.h"
#include "lensfitnessobject.h"
#include <iostream>
#include <limits>
#include <assert.h>

#include "debugnew.h"

using namespace std;

namespace grale
{

GridLensInversionGenomeBase::GridLensInversionGenomeBase(LensInversionGAFactoryCommon *f, int numBasisFunctions, int numSheets)
{
	m_pFactory = f;
	m_fitnessComp = 0;
	
	for (int i = 0 ; i < GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP ; i++)
		m_fitnessValues[i] = 0;
	
	m_scaleFactor = 1.0;
	m_basisFunctionWeights.resize(numBasisFunctions);

	SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());

	for (int i = 0 ; i < numBasisFunctions ; i++)
		m_basisFunctionWeights[i] = (float)uniDist.pickNumber();

	if (m_pFactory->allowNegativeValues())
	{
		for (int i = 0 ; i < numBasisFunctions ; i++)
			m_basisFunctionWeights[i] -= 0.5;
	}

	m_sheetValues.resize(numSheets);
	for (auto &s : m_sheetValues)
		s = uniDist.pickNumber();
}

GridLensInversionGenomeBase::GridLensInversionGenomeBase(LensInversionGAFactoryCommon *f, vector<float> &basisFunctionWeights, vector<float> &sheetValues)
{
	m_pFactory = f;
	m_basisFunctionWeights = std::move(basisFunctionWeights);
	m_sheetValues = std::move(sheetValues);
	
	m_fitnessComp = 0;
	
	for (int i = 0 ; i < GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP ; i++)
		m_fitnessValues[i] = 0;

	m_scaleFactor = 1.0;

#ifndef NDEBUG
	if (!m_pFactory->allowNegativeValues())
	{
		for (auto x : m_basisFunctionWeights)
			assert(x >= 0);
	}
#endif // NDEBUG
}

GridLensInversionGenomeBase::~GridLensInversionGenomeBase()
{
}
	
bool GridLensInversionGenomeBase::calculateFitness()
{
	if (!m_pFactory->calculateFitness(m_basisFunctionWeights, m_sheetValues,
	                                  m_scaleFactor, m_fitnessValues))
		return false;

#ifndef NDEBUG
	if (getenv("GRALE_DEBUG_DUMPFITNESS"))
		m_pFactory->sendMessage("DEBUG: fitness = " + getFitnessDescription());
#endif // !NDEBUG

	return true;
}

mogal::Genome *GridLensInversionGenomeBase::reproduce(const mogal::Genome *g) const
{
	int numBasisFunctions = m_basisFunctionWeights.size();
	std::vector<float> newBasisFunctionsWeights(numBasisFunctions);
	SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());

	const GridLensInversionGenomeBase *pParents[] = { 
		(const GridLensInversionGenomeBase *)g, 
		this 
	};

	auto pickParent = [pParents, &uniDist]()
	{
		float x = (float)uniDist.pickNumber();
		const GridLensInversionGenomeBase *pParent = pParents[(x < 0.5f)?1:0];
		return pParent;
	};

	auto genomeUniformCrossover = [&newBasisFunctionsWeights, numBasisFunctions, pickParent](auto check)
	{
		for (int i = 0 ; i < numBasisFunctions ; i++)
		{
			const GridLensInversionGenomeBase *pParent = pickParent();
			newBasisFunctionsWeights[i] = pParent->m_scaleFactor*pParent->m_basisFunctionWeights[i];
			check(newBasisFunctionsWeights[i]);
		}
	};

	auto noChange = [](float &x) { };
	auto clamp = [](float &x) { if (x < 0) x = 0; };

	if (m_pFactory->allowNegativeValues())
		genomeUniformCrossover(noChange);
	else
		genomeUniformCrossover(clamp);

	vector<float> newSheetValues(m_sheetValues.size());
	for (int i = 0 ; i < newSheetValues.size() ; i++)
		newSheetValues[i] = pickParent()->m_sheetValues[i];

	// This moves newBasisFunctionsWeights and newSheetValues
	return new GridLensInversionGenomeBase(m_pFactory, newBasisFunctionsWeights, newSheetValues);
}

mogal::Genome *GridLensInversionGenomeBase::clone() const
{
	// We need a copy, as the constructor we're going to use will std::move the contents
	vector<float> basisFunctionWeightsCopy { m_basisFunctionWeights };
	vector<float> sheetValuesCopy { m_sheetValues };
	GridLensInversionGenomeBase *g = new GridLensInversionGenomeBase(m_pFactory, basisFunctionWeightsCopy, sheetValuesCopy);

	g->setFitnessValues(m_fitnessValues);
	g->setScaleFactor(m_scaleFactor);

	return g;
}

void GridLensInversionGenomeBase::mutate()
{
	int numBasisFunctions = m_basisFunctionWeights.size();
	float chanceMultiplier = m_pFactory->getChanceMultiplier();
	float chance = chanceMultiplier/((float)numBasisFunctions);
	SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());

	auto absVal = [](auto x) { return ABS(x); };
	auto noChange = [](auto x) { return x; };
	auto getMaxVal = [this, numBasisFunctions](auto fn)
	{
		float maxVal = 0;

		for (int i = 0 ; i < numBasisFunctions ; i++)
		{
			float x = fn(m_basisFunctionWeights[i]);
			if (maxVal < x)
				maxVal = x;
		}
		return maxVal;
	};

	float maxVal = (m_pFactory->allowNegativeValues())?getMaxVal(absVal):getMaxVal(noChange);
	// In the past, this 'maxVal' value was rescaled to 0.5, the unit range then 
	// corresponds to twice this. This means that the unit based mutation amplitude
	// needs to be scaled by 2*maxVal
	float rescale = 2.0f*maxVal;
	float mutationAmplitude = m_pFactory->getMutationAmplitude() * rescale;

	auto chanceSetUniform = [chance, &uniDist, rescale](float &x, float mult, float offset)
	{
		if ((float)uniDist.pickNumber() < chance)
			x = ((float)uniDist.pickNumber()*mult - offset)*rescale;
	};

	auto chanceSetUniformAllScalableBasisFunctions = [this, chance, &uniDist, numBasisFunctions, chanceSetUniform](float mult, float offset)
	{
		for (int i = 0 ; i < numBasisFunctions ; i++)
			chanceSetUniform(m_basisFunctionWeights[i], mult, offset);
	};

	auto chanceSetSmallDiff = [chance, &uniDist, mutationAmplitude](float &target, float yMin, float yMax)
	{
		if ((float)uniDist.pickNumber() < chance)
		{
			// allow larger mutations with smaller probablility
			// p(x) = (2/Pi)*1/(x^2+1)
			// cfr anomalous diffusion
			float p = (float)uniDist.pickNumber()*2.0f-1.0f;
			float x = TAN(p*(float)(CONST_PI/4.0))*mutationAmplitude;
			float y = x+target;

			if (y < yMin)
				y = yMin;
			else if (y > yMax)
				y = yMax;
			
			target = y;
		}
	};

	auto chanceSetSmallDiffAllScalableBasisFunctions = [this, chance, &uniDist, numBasisFunctions, chanceSetSmallDiff](float yMin, float yMax)
	{
		for (int i = 0 ; i < numBasisFunctions ; i++)
			chanceSetSmallDiff(m_basisFunctionWeights[i], yMin, yMax);
	};

	if (m_pFactory->useAbsoluteMutation())
	{
		if (m_pFactory->allowNegativeValues())
			chanceSetUniformAllScalableBasisFunctions(2.0f, 1.0f);
		else
			chanceSetUniformAllScalableBasisFunctions(1.0f, 0.0f);
	}
	else
	{
		if (m_pFactory->allowNegativeValues())
			chanceSetSmallDiffAllScalableBasisFunctions(-rescale, rescale);
		else
			chanceSetSmallDiffAllScalableBasisFunctions(0.0f, rescale);
	}

	if (m_pFactory->useAbsoluteMutation())
	{
		for (auto &v : m_sheetValues)
			chanceSetUniform(v, 1.0f, 0.0f);
	}
	else
	{
		for (auto &v : m_sheetValues)
			chanceSetSmallDiff(v, 0.0f, 1.0f);
	}
}

std::string GridLensInversionGenomeBase::getFitnessDescription() const
{
	char str[1024];
	char *p = str;
	int num = getFactory()->getNumberOfFitnessComponents();

	for (int i = 0 ; i < num ; i++)
	{
		sprintf(p,"%g ",(double)m_fitnessValues[i]);
		p += strlen(p);
	}
	return std::string(str);
}

} // end namespace
