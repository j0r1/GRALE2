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

#ifndef GRALE_LENSFITNESSGENERAL_H

#define GRALE_LENSFITNESSGENERAL_H

#include <grale/lensfitnessobject.h>

namespace grale
{

class FitnessComponent;
class FitnessComponentCache;

class LensFitnessGeneral : public LensFitnessObject
{
public:
	LensFitnessGeneral();
	~LensFitnessGeneral();
	bool init(double z_d, std::list<ImagesDataExtended *> &images, 
	          std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams) override;
	std::string getUsage() const override;

	ConfigurationParameters *getDefaultParametersInstance() const override;
	std::string getFitnessComponentsDescription() const						{ return m_fitnessComponentDescription; }
	int getImagesGroupSize() const											{ return 1; }

	bool shortNeedInverseMagnifications() const 							{ return m_shortInverse; }
	bool shortNeedShearComponents() const									{ return m_shortShear; }
	bool shortNeedConvergence() const										{ return m_shortConvergence; }
	const std::vector<bool> *getShortInverseMagnificationFlags() const		{ return &m_shortInverseFlags; }
	const std::vector<bool> *getShortShearComponentFlags() const 			{ return &m_shortShearFlags; }
	const std::vector<bool> *getShortConvergenceFlags() const				{ return &m_shortConvergenceFlags; }

	bool totalNeedInverseMagnifications() const								{ return m_totalInverse; }
	bool totalNeedShearComponents() const									{ return m_totalShear; }
	bool totalNeedConvergence() const										{ return m_totalConvergence; }
	const std::vector<bool> *getTotalInverseMagnificationFlags() const		{ return &m_totalInverseFlags; }
	const std::vector<bool> *getTotalShearComponentFlags() const			{ return &m_totalShearFlags; }
	const std::vector<bool> *getTotalConvergenceFlags() const 				{ return &m_totalConvergenceFlags; }

	void getTotalCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
			       std::vector<bool> &potentialFlags) const					{ deflectionFlags = m_totalDeflectionFlags; derivativeFlags = m_totalDerivativeFlags; potentialFlags = m_totalPotentialFlags; }
	void getTotalStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const 
																			{ *pStoreIntens = m_totalStoreIntens; *pStoreTimeDelay = m_totalStoreTimeDelay; *pStoreShearInfo = m_totalStoreShear; }

	void getShortCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
			       std::vector<bool> &potentialFlags) const					{ deflectionFlags = m_shortDeflectionFlags; derivativeFlags = m_shortDerivativeFlags; potentialFlags = m_shortPotentialFlags; }
	void getShortStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const	
																			{ *pStoreIntens = m_shortStoreIntens; *pStoreTimeDelay = m_shortStoreTimeDelay; *pStoreShearInfo = m_shortStoreShear; }

	int getNumberOfFitnessComponents() const								{ return m_numFitnessComponents; }

	bool calculateMassScaleFitness(const ProjectedImagesInterface &inf, float &fitness) const override;
	bool calculateOverallFitness(const ProjectedImagesInterface &inf, float *fitnessvalues) const override;
private:
	void clear();

	bool m_initialized;
	int m_numFitnessComponents;

	bool m_totalInverse, m_totalShear, m_totalConvergence;
	bool m_totalStoreIntens, m_totalStoreTimeDelay, m_totalStoreShear;
	bool m_shortInverse, m_shortShear, m_shortConvergence;
	bool m_shortStoreIntens, m_shortStoreTimeDelay, m_shortStoreShear;
	std::vector<bool> m_totalDeflectionFlags, m_totalDerivativeFlags, m_totalPotentialFlags;
	std::vector<bool> m_totalInverseFlags, m_totalShearFlags, m_totalConvergenceFlags;
	std::vector<bool> m_shortDeflectionFlags, m_shortDerivativeFlags, m_shortPotentialFlags;
	std::vector<bool> m_shortInverseFlags, m_shortShearFlags, m_shortConvergenceFlags;

	std::vector<FitnessComponent *> m_totalComponents, m_calculationOrderComponents;
	FitnessComponent *m_pShortComponent;
	bool m_deleteShortComponent;

	FitnessComponentCache *m_pCache;
	std::string m_fitnessComponentDescription;
};

} // end namespace

#endif // GRALE_LENSFITNESSGENERAL_H

