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

#ifndef GRALE_GRIDLENSINVERSIONGENOMEBASE_H

#define GRALE_GRIDLENSINVERSIONGENOMEBASE_H

#include "graleconfig.h"
#include "mathfunctions.h"
#include <mogal/genome.h>
#include <string.h>
#include <string>
#include <vector>

#define GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP 16

namespace grale
{
	
class LensInversionGAFactoryCommon;
class GravitationalLens;	

class GRALE_IMPORTEXPORT GridLensInversionGenomeBase : public mogal::Genome
{
public:
	GridLensInversionGenomeBase(LensInversionGAFactoryCommon *f, int nummasses, int numSheets);
	// This will move the contents of masses and sheetValues!
	GridLensInversionGenomeBase(LensInversionGAFactoryCommon *f, std::vector<float> &masses, std::vector<float> &sheetValues);
	~GridLensInversionGenomeBase();
	
	void setActiveFitnessComponent(int i)								{ m_fitnessComp = i; }
	bool calculateFitness();
	bool isFitterThan(const mogal::Genome *g) const							{ const GridLensInversionGenomeBase *g2 = (const GridLensInversionGenomeBase *)g; if (m_fitnessValues[m_fitnessComp] < g2->getFitnessValues()[m_fitnessComp]) return true; return false; }
	mogal::Genome *reproduce(const mogal::Genome *g) const;
	Genome *clone() const;
	void mutate();
	std::string getFitnessDescription() const;
		
	const float *getFitnessValues() const								{ return m_fitnessValues; }
	float getScaleFactor() const									{ return m_scaleFactor; }
	void setFitnessValues(const float *f)								{ for (int i = 0 ; i < GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP ; i++) m_fitnessValues[i] = f[i]; }
	void setScaleFactor(float s)									{ m_scaleFactor = s; }
	
	const std::vector<float> &getMasses() const							{ return m_masses; }
	const std::vector<float> &getSheetValues() const					{ return m_sheetValues; }

	void copyFitnessValuesTo(float *pDestination, int amount) const					{ memcpy(pDestination, m_fitnessValues, sizeof(float)*amount); }
protected:
	LensInversionGAFactoryCommon *getFactory() const						{ return m_pFactory; }
private:
	int m_fitnessComp;
	float m_fitnessValues[GRIDLENSINVERSIONGENOMEBASE_MAXFITNESSCOMP];	
	float m_scaleFactor;
	std::vector<float> m_masses;
	std::vector<float> m_sheetValues;

	LensInversionGAFactoryCommon *m_pFactory;
};

} // end namespace

#endif // GRALE_GRIDLENSINVERSIONGENOMEBASE_H

