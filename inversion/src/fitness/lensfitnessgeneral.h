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

#pragma once

#include "graleconfig.h"
#include "lensfitnessobject.h"
#include "fitnesscomponent.h"
#include "cosmology.h"
#include <memory>
#include <set>

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

	std::unique_ptr<ConfigurationParameters> getDefaultParametersInstance() const override;
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

	void getTotalCalcFlags(std::vector<bool> &deflectionFlags, 
	                       std::vector<bool> &derivativeFlags, 
	                       std::vector<bool> &potentialFlags,
						   std::vector<bool> &secondDerivFlags) const override { deflectionFlags = m_totalDeflectionFlags; derivativeFlags = m_totalDerivativeFlags; potentialFlags = m_totalPotentialFlags; secondDerivFlags = m_totalSecondDerivFlags; }

	void getShortCalcFlags(std::vector<bool> &deflectionFlags,
	                       std::vector<bool> &derivativeFlags, 
	                       std::vector<bool> &potentialFlags,
						   std::vector<bool> &secondDerivFlags) const override { deflectionFlags = m_shortDeflectionFlags; derivativeFlags = m_shortDerivativeFlags; potentialFlags = m_shortPotentialFlags; secondDerivFlags = m_shortSecondDerivFlags; }

	int getNumberOfFitnessComponents() const								{ return m_numFitnessComponents; }

	bool isNegativeLogProb_Short() const override { assert(m_pShortComponent); return m_pShortComponent->isNegativeLogProb(); }
	bool isNegativeLogProb_Overall(int comp) const override { assert(comp >= 0 && comp < m_totalComponents.size()); assert(m_totalComponents[comp]); return m_totalComponents[comp]->isNegativeLogProb(); }

	bool calculateMassScaleFitness(const ProjectedImagesInterface &inf, float &fitness) const override;
	bool calculateOverallFitness(const ProjectedImagesInterface &inf, float *fitnessvalues) const override;
private:
	void clear();
	bool processGeneralParameters(const ConfigurationParameters *pParams);
	bool processComponentParameters(const ConfigurationParameters *pParams);
	bool setFitnessOptions(FitnessComponent &pComp, const ConfigurationParameters *pParams);
	std::set<std::string> getSupportedTypeNames();
	bool checkImagesDataParameters(std::list<ImagesDataExtended *> &images);
	bool inspectImagesByComponents(
		std::list<ImagesDataExtended *> &images,
		std::vector<std::shared_ptr<FitnessComponent>> &components,
		std::vector<bool> &deflectionFlags,
		std::vector<bool> &derivativeFlags,
		std::vector<bool> &potentialFlags,
		std::vector<bool> &secondDerivFlags,
		std::vector<bool> &inverseFlags,
		std::vector<bool> &shearFlags,
		std::vector<bool> &convergenceFlags,
		bool &calcInverse, bool &calcShear, bool &calcConvergence);
	bool checkUnusedImagesDataParameters(std::list<ImagesDataExtended *> &images);
	bool checkAllImagesUsed(std::list<ImagesDataExtended *> &images);
	bool finalizeAndCleanUnusedComponents(float z_d, const Cosmology *pCosmology);
	std::unique_ptr<FitnessComponent> totalToShort(const std::vector<std::shared_ptr<FitnessComponent>> &total, std::vector<int> &shortImageIndices,
							       const ConfigurationParameters &params);
	void buildShortImagesList(std::list<ImagesDataExtended *> &images,
		const std::vector<int> &shortIndices,
		std::list<ImagesDataExtended *> &shortImages);

	bool m_initialized;
	int m_numFitnessComponents;

	bool m_totalInverse, m_totalShear, m_totalConvergence;
	bool m_shortInverse, m_shortShear, m_shortConvergence;
	std::vector<bool> m_totalDeflectionFlags, m_totalDerivativeFlags, m_totalPotentialFlags, m_totalSecondDerivFlags;
	std::vector<bool> m_totalInverseFlags, m_totalShearFlags, m_totalConvergenceFlags;
	std::vector<bool> m_shortDeflectionFlags, m_shortDerivativeFlags, m_shortPotentialFlags, m_shortSecondDerivFlags;
	std::vector<bool> m_shortInverseFlags, m_shortShearFlags, m_shortConvergenceFlags;

	std::vector<std::shared_ptr<FitnessComponent>> m_totalComponents, m_calculationOrderComponents;
	std::shared_ptr<FitnessComponent> m_pShortComponent;

	std::shared_ptr<FitnessComponentCache> m_pCache;
	std::string m_fitnessComponentDescription;

	std::shared_ptr<Cosmology> m_cosmology;
};

} // end namespace
