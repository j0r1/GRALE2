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
#include "gridlensinversiongafactorybase.h"
#include "gridlensinversiongafactoryparams.h"
#include "multipleplummerlens.h"
#include "multiplesquarelens.h"
#include "multiplegausslens.h"
#include "simpleuniformdistribution.h"
#include "constants.h"
#include "lensfitnessobject.h"
#include "deflectionmatrix.h"
#include "backprojectmatrixnew.h"
#include <iostream>
#include <limits>

#include "debugnew.h"

using namespace std;

namespace grale
{

GridLensInversionGenomeBase::GridLensInversionGenomeBase(GridLensInversionGAFactoryBase *f, int nummasses, bool useGenomeMassSheet)
{
	m_pFactory = f;
	m_fitnessComp = 0;
	m_sheetFactor = 0;
	
	for (int i = 0 ; i < GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP ; i++)
		m_fitnessValues[i] = 0;
	
	m_scaleFactor = 1.0;
	m_masses.resize(nummasses);

	SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());

	for (int i = 0 ; i < nummasses ; i++)
		m_masses[i] = (float)uniDist.pickNumber();

	if (m_pFactory->allowNegativeValues())
	{
		for (int i = 0 ; i < nummasses ; i++)
			m_masses[i] -= 0.5;
	}

	if (useGenomeMassSheet)
		m_sheetValue = uniDist.pickNumber();
	else
		m_sheetValue = -1;

	rescale();
}

GridLensInversionGenomeBase::GridLensInversionGenomeBase(GridLensInversionGAFactoryBase *f, const std::vector<float> &masses, float sheetValue)
{
	m_pFactory = f;
	m_masses = masses;
	m_sheetValue = sheetValue;
	m_sheetFactor = 0;
	
	m_fitnessComp = 0;
	
	for (int i = 0 ; i < GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP ; i++)
		m_fitnessValues[i] = 0;

	m_scaleFactor = 1.0;

	rescale();
}

GridLensInversionGenomeBase::~GridLensInversionGenomeBase()
{
}
	
static float LogTrans(float x) { return LN(x); }
static float ExpTrans(float x) { return EXP(x); }
static float IdentityTrans(float x) { return x; }

bool GridLensInversionGenomeBase::calculateFitness()
{
	int nummasses = m_masses.size();
	float sheetFactor = 0;
	float startValue, startValue0;
	float stopValue, stopValue0;
	float sheetStartValue, sheetStartValue0;
	float sheetStopValue, sheetStopValue0;

	int numiterations;
	int numiterationsteps;
	int numiterationsteps2;

	initializeNewCalculation();

	// TODO: adjust mass scale depending on sheet value?

	if (m_pFactory->allowNegativeValues())
	{
		// TODO: is this enough? or do we need some additional constraints?

		float massum = 0;
		const float *weights = m_pFactory->getMassWeights();
		
		// TODO: should find something better!!
		for (int i = 0 ; i < nummasses ; i++)
			massum += ABS(m_masses[i]*weights[i]);

		m_pFactory->getGenomeCalculationParameters(startValue, stopValue, numiterationsteps, numiterations, numiterationsteps2);
		m_pFactory->getGenomeSheetCalculationParameters(sheetStartValue, sheetStopValue);

		// we want a scale factor of 1 to correspond to the total mass specified by the
		// massscale used in the BackProjectMatrix

		startValue /= massum;
		stopValue /= massum;
	}
	else
	{
		float massum = 0;
		const float *weights = m_pFactory->getMassWeights();

		for (int i = 0 ; i < nummasses ; i++)
			massum += m_masses[i]*weights[i];

		m_pFactory->getGenomeCalculationParameters(startValue, stopValue, numiterationsteps, numiterations, numiterationsteps2);
		m_pFactory->getGenomeSheetCalculationParameters(sheetStartValue, sheetStopValue);

		// we want a scale factor of 1 to correspond to the total mass specified by the
		// massscale used in the BackProjectMatrix

		startValue /= massum;
		stopValue /= massum;
	}

	// Forward transform
	float (*FT)(float x) = IdentityTrans;
	// Inverse transform
	float (*IT)(float x) = IdentityTrans;

	if (m_pFactory->useLogarithmicScaleSearch())
	{
		FT = LogTrans;
		IT = ExpTrans;
	}

	startValue = FT(startValue);
	stopValue = FT(stopValue);

	// Store these start values so that we don't do out of bounds
	startValue0 = startValue;
	stopValue0 = stopValue;

	if (1) // (!loopSheet)
	{
		float currentBestFitness = std::numeric_limits<float>::max();
		float currentBestScaleFactor = 1.0f;

		sheetFactor = m_pFactory->getSheetScale()*m_sheetValue;		

 		for (int i = 0 ; i < numiterations ; i++)
		{
			float stepsize = (stopValue-startValue)/((float)(numiterationsteps-1));
			float s = startValue;

			for (int j = 0 ; j < numiterationsteps ; j++, s += stepsize)
			{
				float realScale = IT(s);
				float f;
				
				if (!calculateMassScaleFitness(realScale, sheetFactor, f))
					return false; // error should already have been set

				if (f < currentBestFitness)
				{
					currentBestFitness = f;
					currentBestScaleFactor = s;
				}
			}

			startValue = currentBestScaleFactor-stepsize;
			stopValue = currentBestScaleFactor+stepsize;

			// Make sure we stay within bounds
			if (startValue < startValue0)
				startValue = startValue0;
			if (stopValue > stopValue0)
				stopValue = stopValue0;

			// After the first loop, we're going to just zoom in on the first located value
			// This can use a different number of steps
			numiterationsteps = numiterationsteps2;
		}

		m_scaleFactor = IT(currentBestScaleFactor);
	}

	// Do the final (possibly multi-component) fitness evaluation
	calculateTotalFitness(m_scaleFactor, sheetFactor, m_fitnessValues);

	return true;
}

mogal::Genome *GridLensInversionGenomeBase::reproduce(const mogal::Genome *g) const
{
	int nummasses = m_masses.size();
	std::vector<float> newmasses(nummasses);
	const GridLensInversionGenomeBase *g2 = (const GridLensInversionGenomeBase *)g;
	SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());

	if (m_pFactory->allowNegativeValues())
	{
		for (int i = 0 ; i < nummasses ; i++)
		{
			float x = (float)uniDist.pickNumber();

			if (x < 0.5f)
				newmasses[i] = m_scaleFactor*m_masses[i];
			else
				newmasses[i] = (g2->m_scaleFactor)*(g2->m_masses[i]);
		}
	}
	else
	{
		for (int i = 0 ; i < nummasses ; i++)
		{
			float x = (float)uniDist.pickNumber();

			if (x < 0.5f)
				newmasses[i] = m_scaleFactor*m_masses[i];
			else
				newmasses[i] = (g2->m_scaleFactor)*(g2->m_masses[i]);
			if (newmasses[i] < 0)
				newmasses[i] = 0;
		}
	}

	float newSheetValue = -1;

	if (m_sheetValue >= 0)
	{
		float x = (float)uniDist.pickNumber();

		if (x < 0.5f)
			newSheetValue = m_sheetValue;
		else
			newSheetValue = g2->m_sheetValue;
	}

	return new GridLensInversionGenomeBase(m_pFactory, newmasses, newSheetValue);
}

mogal::Genome *GridLensInversionGenomeBase::clone() const
{
	GridLensInversionGenomeBase *g = new GridLensInversionGenomeBase(m_pFactory, m_masses, m_sheetValue);

	g->setFitnessValues(m_fitnessValues);
	g->setScaleFactor(m_scaleFactor);
	g->setSheetFactor(m_sheetFactor);

	return g;
}

void GridLensInversionGenomeBase::mutate()
{
	int nummasses = m_masses.size();

	if (m_pFactory->allowNegativeValues())
	{
		float chanceMultiplier = m_pFactory->getChanceMultiplier();
		float chance = chanceMultiplier/((float)nummasses);
		SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());
		
		if (m_pFactory->useAbsoluteMutation())
		{
			for (int i = 0 ; i < nummasses ; i++)
			{
				if ((float)uniDist.pickNumber() < chance)
					m_masses[i] = (float)uniDist.pickNumber()*2.0f-1.0f;
			}
		}
		else
		{
			float mutationAmplitude = m_pFactory->getMutationAmplitude();

			for (int i = 0 ; i < nummasses ; i++)
			{
				if ((float)uniDist.pickNumber() < chance)
				{
					// allow larger mutations with smaller probablility
					// p(x) = (2/Pi)*1/(x^2+1)
					// cfr anomalous diffusion
					float p = (float)uniDist.pickNumber()*2.0f-1.0f;
					float x = TAN(p*(float)(CONST_PI/4.0))*mutationAmplitude;
					float y = x+m_masses[i];

					if (y < -1.0f)
						y = -1.0f;
					else if (y > 1.0f)
						y = 1.0f;
					
					m_masses[i] = y;
				}
			}
		}
	}
	else
	{
		float chanceMultiplier = m_pFactory->getChanceMultiplier();
		float chance = chanceMultiplier/((float)nummasses);
		SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());
		
		if (m_pFactory->useAbsoluteMutation())
		{
			for (int i = 0 ; i < nummasses ; i++)
			{
				if ((float)uniDist.pickNumber() < chance)
					m_masses[i] = (float)uniDist.pickNumber();
			}
		}
		else
		{
			float mutationAmplitude = m_pFactory->getMutationAmplitude();

			for (int i = 0 ; i < nummasses ; i++)
			{
				if ((float)uniDist.pickNumber() < chance)
				{
					// allow larger mutations with smaller probablility
					// p(x) = (2/Pi)*1/(x^2+1)
					// cfr anomalous diffusion
					float p = (float)uniDist.pickNumber()*2.0f-1.0f;
					float x = TAN(p*(float)(CONST_PI/4.0))*mutationAmplitude;
					float y = x+m_masses[i];

					if (y < 0)
						y = 0;
					else if (y > 1.0)
						y = 1.0;
					
					m_masses[i] = y;
				}
			}
		}
	}

	if (m_sheetValue >= 0)
	{
		float chanceMultiplier = m_pFactory->getChanceMultiplier();
		float chance = chanceMultiplier/((float)nummasses);
		SimpleUniformDistribution uniDist(m_pFactory->getRandomNumberGenerator());
		
		if (m_pFactory->useAbsoluteMutation())
		{
			if ((float)uniDist.pickNumber() < chance)
				m_sheetValue = (float)uniDist.pickNumber();
		}
		else
		{
			float mutationAmplitude = m_pFactory->getMutationAmplitude();

			if ((float)uniDist.pickNumber() < chance)
			{
				// allow larger mutations with smaller probablility
				// p(x) = (2/Pi)*1/(x^2+1)
				// cfr anomalous diffusion
				float p = (float)uniDist.pickNumber();
				float x = TAN(p*(float)(CONST_PI/4.0))*mutationAmplitude;
				float y = x+m_sheetValue;

				if (y < 0)
					y = 0;
				else if (y > 1.0f)
					y = 1.0f;
				
				m_sheetValue = y;
			}
		}
	}

	rescale();
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
	
void GridLensInversionGenomeBase::rescale()
{
	int nummasses = m_masses.size();

	if (m_pFactory->allowNegativeValues())
	{
		float maxval = 0;
		
		for (int i = 0 ; i < nummasses ; i++)
		{
			if (ABS(maxval) < ABS(m_masses[i]))
				maxval = ABS(m_masses[i]);
		}

		if (maxval > 0)
		{
			float f = 0.5f/maxval;
			for (int i = 0 ; i < nummasses ; i++)
				m_masses[i] *= f;
		}
	}
	else
	{
		float maxval = 0;
		for (int i = 0 ; i < nummasses ; i++)
			maxval = MAX(maxval, m_masses[i]);
		
		if (maxval > 0)
		{
			float f = 0.5f/maxval;
			for (int i = 0 ; i < nummasses ; i++)
			{
				m_masses[i] *= f;
				if (m_masses[i] < 0)
					m_masses[i] = 0;
			}
		}
	}
}

GravitationalLens *GridLensInversionGenomeBase::createLens(double *totalmass, std::string &errstr) const
{
	return m_pFactory->createLens(m_masses, m_sheetValue, m_scaleFactor, totalmass, errstr);
}

void GridLensInversionGenomeBase::initializeNewCalculation()
{
	GridLensInversionGAFactoryBase *pFactory = getFactory();
	DeflectionMatrix *pDeflectionMatrix = pFactory->getDeflectionMatrix();
	BackProjectMatrixNew *pTotalBackProjectMatrix = pFactory->getTotalBackProjectMatrix();
	BackProjectMatrixNew *pShortBackProjectMatrix = pFactory->getShortBackProjectMatrix();

	pDeflectionMatrix->calculateBasisMatrixProducts(getMasses(), true, true, true);
	pTotalBackProjectMatrix->storeDeflectionMatrixResults();
	pShortBackProjectMatrix->storeDeflectionMatrixResults();
}

bool GridLensInversionGenomeBase::calculateMassScaleFitness(float scaleFactor, float sheetScale, float &fitness)
{
	GridLensInversionGAFactoryBase *pFactory = getFactory();
	BackProjectMatrixNew *pShortBackProjectMatrix = pFactory->getShortBackProjectMatrix();
	LensFitnessObject &fitnessFunction = *(pFactory->getFitnessObject());

	pShortBackProjectMatrix->calculate(scaleFactor, sheetScale);
	if (fitnessFunction.shortNeedInverseMagnifications())
		pShortBackProjectMatrix->calculateInverseMagnifications(*(fitnessFunction.getShortInverseMagnificationFlags()));
	if (fitnessFunction.shortNeedShearComponents())
		pShortBackProjectMatrix->calculateShearComponents(*(fitnessFunction.getShortShearComponentFlags()));
	if (fitnessFunction.shortNeedConvergence())
		pShortBackProjectMatrix->calculateConvergence(*(fitnessFunction.getShortConvergenceFlags()));

	if (!fitnessFunction.calculateMassScaleFitness(*pShortBackProjectMatrix, fitness))
	{
		cerr << "ERROR: Unable to calculate mass scale fitness: " << fitnessFunction.getErrorString() << endl;
		return false;
	}
	return true;
}

bool GridLensInversionGenomeBase::calculateTotalFitness(float scaleFactor, float sheetScale, float *pFitnessValues)
{
	GridLensInversionGAFactoryBase *pFactory = getFactory();
	BackProjectMatrixNew *pTotalBackProjectMatrix = pFactory->getTotalBackProjectMatrix();
	LensFitnessObject &fitnessFunction = *(pFactory->getFitnessObject());
	
	pTotalBackProjectMatrix->calculate(scaleFactor, sheetScale);
	if (fitnessFunction.totalNeedInverseMagnifications())
		pTotalBackProjectMatrix->calculateInverseMagnifications(*(fitnessFunction.getTotalInverseMagnificationFlags()));
	if (fitnessFunction.totalNeedShearComponents())
		pTotalBackProjectMatrix->calculateShearComponents(*(fitnessFunction.getTotalShearComponentFlags()));
	if (fitnessFunction.totalNeedConvergence())
		pTotalBackProjectMatrix->calculateConvergence(*(fitnessFunction.getTotalConvergenceFlags()));

	if (!fitnessFunction.calculateOverallFitness(*pTotalBackProjectMatrix, pFitnessValues))
	{
		cerr << "ERROR: Unable to calculate full fitness: " << fitnessFunction.getErrorString() << endl;
		return false;
	}
	return true;
}

} // end namespace


