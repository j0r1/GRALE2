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

#ifndef GRALE_LENSFITNESSJ1004_H

#define GRALE_LENSFITNESSJ1004_H

#include <grale/lensfitnessobject.h>
#include <grale/rectangle2d.h>
#include <grale/triangleindices.h>
#include "pointgroupstorage.h"
#include <vector>
#include <list>

namespace grale
{

class FitnessComponentCache;

class LensFitnessJ1004Test : public LensFitnessObject
{
public:
	LensFitnessJ1004Test();
	~LensFitnessJ1004Test();
	bool init(double z_d, std::list<ImagesDataExtended *> &images, 
	          std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams) override;
	std::string getUsage() const override;

	bool isInitialized() const									{ return m_initialized; }
	int getImagesGroupSize() const									{ return 3; }
	int getNumberOfFitnessComponents() const							{ return 4; }

	bool shortNeedInverseMagnifications() const 							{ return false; }
	bool shortNeedShearComponents() const								{ return false; }
	bool shortNeedConvergence() const 								{ return false; }
	const std::vector<bool> *getShortInverseMagnificationFlags() const				{ return 0; }
	const std::vector<bool> *getShortShearComponentFlags() const 					{ return 0; }
	const std::vector<bool> *getShortConvergenceFlags() const					{ return 0; }
	bool totalNeedInverseMagnifications() const							{ return true; }
	bool totalNeedShearComponents() const 								{ return false; }
	bool totalNeedConvergence() const 								{ return false; }
	const std::vector<bool> *getTotalInverseMagnificationFlags() const				{ return &m_totalInverseMagnificationFlags; }
	const std::vector<bool> *getTotalShearComponentFlags() const					{ return 0; }
	const std::vector<bool> *getTotalConvergenceFlags() const 					{ return 0; }

	void getTotalCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
			       std::vector<bool> &potentialFlags) const					{ deflectionFlags = m_totalDeflectionFlags; derivativeFlags = m_totalDerivativeFlags; potentialFlags = m_totalPotentialFlags; }
	void getTotalStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const { *pStoreIntens = false; *pStoreTimeDelay = true; *pStoreShearInfo = false; }
	void getShortCalcFlags(std::vector<bool> &deflectionFlags, std::vector<bool> &derivativeFlags, 
			       std::vector<bool> &potentialFlags) const					{ deflectionFlags = m_shortDeflectionFlags; derivativeFlags = m_shortDerivativeFlags; potentialFlags = m_shortPotentialFlags; }
	void getShortStoreFlags(bool *pStoreIntens, bool *pStoreTimeDelay, bool *pStoreShearInfo) const	{ *pStoreIntens = false; *pStoreTimeDelay = false; *pStoreShearInfo = false; }

	bool calculateMassScaleFitness(const ProjectedImagesInterface &inf, float &fitness) const override;
	bool calculateOverallFitness(const ProjectedImagesInterface &inf, float *fitnessvalues) const override;
private:
	int storeTriangleLine(int sourcePos, int imgIndex, int point1, int point2);

	std::vector<std::vector<TriangleIndices> > m_nullTriangles;
	std::vector<std::vector<double> > m_nullTriangleAreas;
	std::vector<std::vector<std::vector<std::pair<int,int> > > > m_lineSegments;
	mutable std::vector<std::vector<std::vector<bool> > > m_lineSegmentFlags;
	mutable std::vector<std::vector<std::vector<Vector2D<float> > > > m_lineSegmentIntersections;
	std::vector<std::vector<std::vector<TriangleIndices> > > m_critTriangles;
	std::vector<int> m_fitnessMask;
	int m_halfNumberOfSources;
	bool m_initialized;

	PointGroupStorage m_pointGroupsAll;
	PointGroupStorage m_pointGroupsShort;

	std::vector<bool> m_totalInverseMagnificationFlags;
	std::vector<bool> m_totalDeflectionFlags, m_totalDerivativeFlags, m_totalPotentialFlags;
	std::vector<bool> m_shortDeflectionFlags, m_shortDerivativeFlags, m_shortPotentialFlags;

	std::vector<int> m_sourceIndices, m_nullIndices, m_critSources, m_critIndices, m_shortSourceIndices, m_tdIndices;
	std::vector<bool> m_rectFlags, m_pointGroupFlags;
	std::vector<float> m_nullWeights;
	FitnessComponentCache *m_pCache;
};

} // end namespace

#endif // GRALE_LENSFITNESSJ1004_H

