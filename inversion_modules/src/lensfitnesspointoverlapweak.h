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

#ifndef GRALE_LENSFITNESSPOINTOVERLAPWEAK_H

#define GRALE_LENSFITNESSPOINTOVERLAPWEAK_H

#include <grale/lensfitnessobject.h>
#include <grale/rectangle2d.h>
#include <grale/triangleindices.h>
#include <vector>
#include <list>

namespace grale
{

class LensFitnessPointOverlapNullWeak : public LensFitnessObject
{
public:
	LensFitnessPointOverlapNullWeak();
	~LensFitnessPointOverlapNullWeak();
	
	bool init(double z_d, std::list<ImagesDataExtended *> &images, 
	          std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams) override;
	std::string getUsage() const override;

	bool isInitialized() const																			{ return m_initialized; }

	// Since we're using both strong and weak lensing data, there's no group size
	int getImagesGroupSize() const																		{ return 1; }
	std::string getFitnessComponentsDescription() const	{ return "pointimageoverlap pointimagenull weaklensing"; }

	bool shortNeedInverseMagnifications() const															{ return false; }
	bool shortNeedShearComponents() const																{ return false; }
	bool shortNeedConvergence() const																	{ return false; }
	const std::vector<bool> *getShortInverseMagnificationFlags() const									{ return 0; }
	const std::vector<bool> *getShortShearComponentFlags() const										{ return 0; }
	const std::vector<bool> *getShortConvergenceFlags() const											{ return 0; }

	bool totalNeedInverseMagnifications() const															{ return false; }
	bool totalNeedShearComponents() const																{ return true; }
	bool totalNeedConvergence() const																	{ return true; }
	const std::vector<bool> *getTotalInverseMagnificationFlags() const									{ return 0; }
	const std::vector<bool> *getTotalShearComponentFlags() const										{ return &m_weakOneFlags; }
	const std::vector<bool> *getTotalConvergenceFlags() const											{ return &m_weakOneFlags; }

	void getTotalCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
	                       std::vector<bool> &potentialFlags) const										{ deflectionFlags = m_deflectionFlags; derivativeFlags = m_derivativeFlags; potentialFlags = m_potentialFlags; }
	void getTotalStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const		{ *pStoreIntens = false; *pStoreTimeDelay = false; *pStoreShearInfo = true; }
	void getShortCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
	                       std::vector<bool> &potentialFlags) const										{ deflectionFlags = m_shortDeflectionFlags; derivativeFlags = m_shortDerivativeFlags; potentialFlags = m_shortPotentialFlags; }
	void getShortStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const		{ *pStoreIntens = false; *pStoreTimeDelay = false; *pStoreShearInfo = false; }

	int getNumberOfFitnessComponents() const															{ return 3; }

	bool calculateMassScaleFitness(const ProjectedImagesInterface &interface, float &fitness) const override;
	bool calculateOverallFitness(const ProjectedImagesInterface &interface, float *pFitnessValues) const override;
private:
	bool m_initialized;

	std::vector<std::vector<TriangleIndices> > m_nullTriangles;
	std::vector<bool> m_deflectionFlags, m_derivativeFlags, m_potentialFlags;
	std::vector<bool> m_shortDeflectionFlags, m_shortDerivativeFlags, m_shortPotentialFlags;
	std::vector<float> m_distanceFractions;

	std::vector<bool> m_weakOneFlags;
	int m_weakStartPos;

	std::vector<int> m_sourceIndices, m_shortSourceIndices;
	std::vector<int> m_nullIndices;
	std::vector<float> m_nullWeights;

	std::vector<int> m_weakIndices;
	std::vector<bool> m_reducedShear;
	std::vector<double> m_thresholds;
};

} // end namespace

#endif // GRALE_LENSFITNESSPOINTOVERLAPWEAK_H

