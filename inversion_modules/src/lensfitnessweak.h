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

#ifndef GRALE_LENSFITNESSWEAK_H

#define GRALE_LENSFITNESSWEAK_H

#include "fitnessutil.h"
#include <grale/lensfitnessobject.h>

namespace grale
{

class LensFitnessWeak : public LensFitnessObject
{
public:
	LensFitnessWeak();
	~LensFitnessWeak();

	bool init(double z_d, std::list<ImagesDataExtended *> &images, 
	          std::list<ImagesDataExtended *> &massScaleImages, const ConfigurationParameters *pParams) override;
	std::string getUsage() const override;
	std::string getFitnessComponentsDescription() const													{ return "weaklensing"; }

	void getTotalCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
	                       std::vector<bool> &potentialFlags) const										{ deflectionFlags = m_zeroVector; derivativeFlags = m_oneVector; potentialFlags = m_zeroVector; }
	void getTotalStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const		{ *pStoreIntens = false; *pStoreTimeDelay = false; *pStoreShearInfo = true; }
	 
	int getNumberOfFitnessComponents() const															{ return 1; }
	int getImagesGroupSize () const																		{ return 1; }

	bool totalNeedInverseMagnifications() const															{ return false; }
	bool totalNeedShearComponents() const																{ return true; }
	bool totalNeedConvergence() const																	{ return true; } // needed for reduced shear
	const std::vector<bool> *getTotalInverseMagnificationFlags() const									{ return 0; }
	const std::vector<bool> *getTotalShearComponentFlags() const										{ return &m_oneVector; }
	const std::vector<bool> *getTotalConvergenceFlags() const											{ return &m_oneVector; }

	void getShortCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
	                       std::vector<bool> &potentialFlags) const										{ deflectionFlags = m_zeroVector; derivativeFlags = m_oneVector; potentialFlags = m_zeroVector; }
	void getShortStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const		{ *pStoreIntens = false; *pStoreTimeDelay = false; *pStoreShearInfo = true; }
	bool shortNeedInverseMagnifications() const															{ return false; }
	bool shortNeedShearComponents() const																{ return true; }
	bool shortNeedConvergence() const																	{ return true; } // needed for reduced shear
	const std::vector<bool> *getShortInverseMagnificationFlags() const									{ return 0; }
	const std::vector<bool> *getShortShearComponentFlags() const										{ return &m_oneVector; }
	const std::vector<bool> *getShortConvergenceFlags() const											{ return &m_oneVector; }


	bool calculateMassScaleFitness(const ProjectedImagesInterface &interface, float &fitness) const override { float value[1]; if (!calculateOverallFitness(interface, value)) return false; fitness = value[0]; return true; }
	bool calculateOverallFitness(const ProjectedImagesInterface &interface, float *pFitnessValues) const override;
private:
	bool m_initialized;

	std::vector<int> m_weakIndices;
	std::vector<bool> m_oneVector, m_zeroVector;
	std::vector<WeakLensingType> m_reducedShear;
	std::vector<double> m_thresholds;
};

} // end namespace

#endif // GRALE_LENSFITNESSWEAK_H

