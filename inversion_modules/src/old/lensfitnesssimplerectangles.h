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

#ifndef GRALE_LENSFITNESSSIMPLERECTANGLES_H

#define GRALE_LENSFITNESSSIMPLERECTANGLES_H

#include <grale/lensfitnessobject.h>
#include <grale/rectangle2d.h>
#include <grale/triangleindices.h>
#include <grale/pointgroupstorage.h>
#include <vector>
#include <list>

namespace grale
{

class LensFitnessSimpleRectangles : public LensFitnessObject
{
public:
	LensFitnessSimpleRectangles() { m_initialized = false; }
	~LensFitnessSimpleRectangles() { }
	bool init(double z_d, std::list<ImagesDataExtended *> &images, 
	          std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams) override;
	std::string getUsage() const override;

	int getImagesGroupSize() const									{ return 1; }
	int getNumberOfFitnessComponents() const							{ return 1; }
	std::string getFitnessComponentsDescription() const				 { return "extendedimageoverlap"; }

	bool shortNeedInverseMagnifications() const 							{ return false; }
	bool shortNeedShearComponents() const								{ return false; }
	bool shortNeedConvergence() const 								{ return false; }
	const std::vector<bool> *getShortInverseMagnificationFlags() const				{ return 0; }
	const std::vector<bool> *getShortShearComponentFlags() const 					{ return 0; }
	const std::vector<bool> *getShortConvergenceFlags() const					{ return 0; }
	bool totalNeedInverseMagnifications() const							{ return false; }
	bool totalNeedShearComponents() const 								{ return false; }
	bool totalNeedConvergence() const 								{ return false; }
	const std::vector<bool> *getTotalInverseMagnificationFlags() const				{ return 0; }
	const std::vector<bool> *getTotalShearComponentFlags() const					{ return 0; }
	const std::vector<bool> *getTotalConvergenceFlags() const 					{ return 0; }

	void getTotalCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
			       std::vector<bool> &potentialFlags) const					{ deflectionFlags = m_deflectionFlags; derivativeFlags = m_derivativeFlags; potentialFlags = m_potentialFlags; }
	void getTotalStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const { *pStoreIntens = false; *pStoreTimeDelay = false; *pStoreShearInfo = false; }
	void getShortCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
			       std::vector<bool> &potentialFlags) const					{ }
	void getShortStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const	{ }


	bool calculateMassScaleFitness(const ProjectedImagesInterface &inf, float &fitness) const override;
	bool calculateOverallFitness(const ProjectedImagesInterface &inf, float *fitnessvalues) const override;
private:
	bool m_initialized;
	std::vector<bool> m_deflectionFlags, m_derivativeFlags, m_potentialFlags;
	PointGroupStorage m_pointGroupInfo;

	std::vector<int> m_sourceIndices;
	std::vector<bool> m_useRectFlags;
	std::vector<bool> m_usePointGroupFlags;
};
	
} // end namespace

#endif // GRALE_LENSFITNESSSIMPLERECTANGLES_H

