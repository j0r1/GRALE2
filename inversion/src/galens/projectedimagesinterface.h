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

/**
 * \file projectedimagesinterface.h
 */

#ifndef GRALE_PROJECTEDIMAGESINTERFACE_H

#define GRALE_PROJECTEDIMAGESINTERFACE_H

#include "graleconfig.h"
#include "vector2d.h"
#include "constants.h"
#include "imagesdata.h"
#include <stdio.h>
#include <assert.h>
#include <vector>

namespace grale
{

class ImagesDataExtended;

/** Base class defining the methods that a LensFitnessObject can use 
 *  to calculate the fitness of a specific lens.
 *  Base class defining the methods that a LensFitnessObject can use 
 *  to calculate the fitness of a specific lens. 
 *
 *  Two implementations
 *  currently exist. The first, and simplest, is the ImagesBackProjector
 *  class. It takes a GravitationalLens object (like a lens model
 *  created with GRALESHELL using for example \c 'lens/new/mplummers')
 *  and a list of images data and calculates a number of properties
 *  (backprojected image points, convergence, shear, ...) and stores
 *  others (measured shear for example). It is this implementation
 *  that is used by the GRALESHELL command \c 'imgdata/fitness'.
 *
 *  The second implementation, BackProjectMatrix, is somewhat more
 *  advanced, but is only meant to be used inside the genetic algorithm.
 *  It is more advanced in the sense that it used information provided
 *  by the LensFitnessObject to make sure that the calculations are done
 *  in an efficient way, only calculating those things that are necessary.
 *  This also means that your program will crash if you try to access 
 *  information that the BackProjectMatrix instance believed was
 *  not necessary to calculate.
 */
class GRALE_IMPORTEXPORT ProjectedImagesInterface
{
protected:
	ProjectedImagesInterface()							{ }
public:
	virtual ~ProjectedImagesInterface()						{ }
	
	/** Returns the number of sources for which data is stored ('source' is
	 *  here just an entry in the images data list, so it could also indicate
	 *  the null space data for a particular object) */
	virtual int getNumberOfSources() const;

	/** Returns the number of images stored for a specific source. */
	virtual int getNumberOfImages(int sourceNumber) const;
	
	/** Returns the total number of points of all images corresponding to a source index. */
	virtual int getNumberOfImagePoints(int sourceNumber) const;

	/** Returns the number of points stored for a particular image of a particular source. */
	virtual int getNumberOfImagePoints(int sourceNumber, int imageNumber) const;

	virtual bool hasOriginalProperty(ImagesData::PropertyName n, int sourceNumber) const;
	virtual const float *getOriginalProperties(ImagesData::PropertyName n, int sourceNumber) const;
	virtual const float *getOriginalProperties(ImagesData::PropertyName n, int sourceNumber, int imageNumber) const;

	/** Returns the number of time delays that were stored for a specific source index. */
	virtual int getOriginalNumberOfTimeDelays(int sourceNumber) const;

	/** For a specific source and specific index (from 0 to getOriginalNumberOfTimeDelays - 1),
	 *  the image index for a time delay value will be stored in \c pImg, the point index in \c pPoint
	 *  and the actual measured time delay in \c pDelay. */
	virtual void getOriginalTimeDelay(int sourceNumber, int index, int *pImg, int *pPoint, float *pDelay) const;

	virtual float getDistanceFraction(int sourcenum) const { return m_distanceFractions[sourcenum]; }





	/** Returns the distance to the lens. */
	virtual double getLensDistance() const = 0;

	/** Returns the redshift of the lens. */
	virtual double getLensRedshift() const = 0;

	/** Returns the angular scale in which calculated backprojected and original point positions are expressed. */
	virtual double getAngularScale() const = 0;

	/** Returns the positions in the source plane for all points of a specific source index. */
	virtual const Vector2D<float> *getBetas(int sourcenum) const = 0;

	/** Returns the positions in the source plane for a specific image of a specific source. */
	virtual const Vector2D<float> *getBetas(int sourcenum, int imagenum) const = 0;

	/** Returns the (scaled) deflection angles all points of a specific source index. */
	virtual const Vector2D<float> *getAlphas(int sourcenum) const = 0;

	/** Returns the (scaled) deflection angles for a specific image of a specific source. */
	virtual const Vector2D<float> *getAlphas(int sourcenum, int imagenum) const = 0;

	/** Returns the (stored) positions in the image plane of all points of a specific source. */
	virtual const Vector2D<float> *getThetas(int sourcenum) const = 0;

	/** Returns the (stored) positions in the image plane of a specific image of a specific source. */
	virtual const Vector2D<float> *getThetas(int sourcenum, int imagenum) const = 0;

	// TODO: make pure virtual?
	virtual bool hasRetracedThetas(int sourceNum) const { return false; }
	virtual const Vector2D<float> *getRetracedThetas(int sourceNum) const { return nullptr; }
	virtual const Vector2D<float> *getRetracedThetas(int sourceNum, int imageNum) const { return nullptr; }

	/** Returns alpha_xx values for all points of a specific source. */
	virtual const float *getDerivativesXX(int sourceNumber) const = 0;

	/** Returns alpha_xx values for the points of a specific image of a specific source. */
	virtual const float *getDerivativesXX(int sourceNumber, int imageNumber) const = 0;

	/** Returns alpha_yy values for all points of a specific source. */
	virtual const float *getDerivativesYY(int sourceNumber) const = 0;

	/** Returns alpha_yy values for the points of a specific image of a specific source. */
	virtual const float *getDerivativesYY(int sourceNumber, int imageNumber) const = 0;

	/** Returns alpha_xy values for all points of a specific source. */
	virtual const float *getDerivativesXY(int sourceNumber) const = 0;

	/** Returns alpha_xy values for the points of a specific image of a specific source. */
	virtual const float *getDerivativesXY(int sourceNumber, int imageNumber) const = 0;

	virtual const float *getSecondDerivativesXXX(int sourceNumber) const = 0;
	virtual const float *getSecondDerivativesXXX(int sourceNumber, int imageNumber) const = 0;
	virtual const float *getSecondDerivativesYYY(int sourceNumber) const = 0;
	virtual const float *getSecondDerivativesYYY(int sourceNumber, int imageNumber) const = 0;
	virtual const float *getSecondDerivativesXXY(int sourceNumber) const = 0;
	virtual const float *getSecondDerivativesXXY(int sourceNumber, int imageNumber) const = 0;
	virtual const float *getSecondDerivativesYYX(int sourceNumber) const = 0;
	virtual const float *getSecondDerivativesYYX(int sourceNumber, int imageNumber) const = 0;


	/** Returns the calculated inverse magnification values for all points of a specific
	 *  source (expressed in termes of the intensity scale) */
	virtual const float *getInverseMagnifications(int sourcenum) const = 0;

	/** Returns the calculated inverse magnification values for a specific image of a
	 *  specific source (expressed in termes of the intensity scale) */
	virtual const float *getInverseMagnifications(int sourcenum, int imagenum) const = 0;

	/** Returns the calculated gamma1 values for all the points of the specified source. */
	virtual const float *getShearComponents1(int sourceNumber) const = 0;

	/** Returns the calculated gamma1 values for a specific image of a specific source. */
	virtual const float *getShearComponents1(int sourceNumber, int imageNumber) const = 0;

	/** Returns the calculated gamma2 values for all the points of the specified source. */
	virtual const float *getShearComponents2(int sourceNumber) const = 0;

	/** Returns the calculated gamma2 values for a specific image of a specific source. */
	virtual const float *getShearComponents2(int sourceNumber, int imageNumber) const = 0;

	/** Returns the calculated convergence (kappa) values for all points of the specified source index. */
	virtual const float *getConvergence(int sourceNumber) const = 0;

	/** Returns the calculated convergence (kappa) values for a specific image of a specific source. */
	virtual const float *getConvergence(int sourceNumber, int imageNumber) const = 0;

	// Lensing potential scale is angular scale squared
	// This way, it makes sense in the time delay formula, where 1/2 (beta-theta)^2 - phi(theta)
	// is used
	//virtual double getLensPotentialScale() const = 0;
	virtual const float *getLensPotential(int sourceNumber) const = 0;
	virtual const float *getLensPotential(int sourceNumber, int imageNumber) const = 0;

	/** For the specified point of an image of a source, calculate the time delay if the backprojected
	 *  position of that point were \c beta. */
	virtual float getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const = 0; // in days!
protected:
	void storeOriginalData(const std::vector<ImagesDataExtended *> &images);

	class TimeDelayPoint
	{
	public:
		TimeDelayPoint(int imageIndex, int pointIndex, float timeDelay)
		{
			m_imageIndex = imageIndex;
			m_pointIndex = pointIndex;
			m_timeDelay = timeDelay;
		}

		int getImageIndex() const						{ return m_imageIndex; }
		int getPointIndex() const						{ return m_pointIndex; }
		float getTimeDelay() const						{ return m_timeDelay; }
	private:
		int m_imageIndex;
		int m_pointIndex;
		float m_timeDelay;
	};

	std::vector<std::vector<int> > m_offsets;
	std::vector<std::vector<int> > m_numPoints;
	std::vector<int> m_numTotalPoints;

	std::vector<float> m_distanceFractions;
	std::vector<std::vector<bool>> m_originalPropertyFlags;
	std::vector<std::vector<std::vector<float>>> m_originalPointProperties;
	std::vector<std::vector<TimeDelayPoint> > m_originalTimeDelayInfo;
};

inline int ProjectedImagesInterface::getNumberOfSources() const
{ 
	return m_offsets.size(); 
}

inline int ProjectedImagesInterface::getNumberOfImages(int sourceNumber) const
{ 
	assert(sourceNumber >= 0 && sourceNumber < m_offsets.size());
	return m_offsets[sourceNumber].size(); 
}
	
inline int ProjectedImagesInterface::getNumberOfImagePoints(int sourceNumber) const
{ 
	assert(sourceNumber >= 0 && sourceNumber < m_numTotalPoints.size());
	return m_numTotalPoints[sourceNumber]; 
}

inline int ProjectedImagesInterface::getNumberOfImagePoints(int sourceNumber, int imageNumber) const
{ 
	assert(sourceNumber >= 0 && sourceNumber < m_numPoints.size());
	assert(imageNumber >= 0 && imageNumber < m_numPoints[sourceNumber].size());
	return m_numPoints[sourceNumber][imageNumber]; 
}

inline bool ProjectedImagesInterface::hasOriginalProperty(ImagesData::PropertyName n, int sourceNumber) const
{
	assert((int)n >= 0 && (int)n < m_originalPropertyFlags.size());
	assert(sourceNumber >= 0 && sourceNumber < m_originalPropertyFlags[(int)n].size());
	
	return m_originalPropertyFlags[(int)n][sourceNumber];
}

inline const float *ProjectedImagesInterface::getOriginalProperties(ImagesData::PropertyName n, int sourceNumber) const
{
	assert((int)n >= 0 && (int)n < m_originalPointProperties.size());
	assert(sourceNumber >= 0 && sourceNumber < m_originalPointProperties[(int)n].size());
	assert(m_originalPointProperties[(int)n][sourceNumber].size() > 0);
	return m_originalPointProperties[(int)n][sourceNumber].data();
}

inline const float *ProjectedImagesInterface::getOriginalProperties(ImagesData::PropertyName n, int sourceNumber, int imageNumber) const
{
	assert((int)n >= 0 && (int)n < m_originalPointProperties.size());
	assert(sourceNumber >= 0 && sourceNumber < m_originalPointProperties[(int)n].size());
	assert(sourceNumber < m_offsets.size());
	assert(imageNumber >= 0 && imageNumber < m_offsets[sourceNumber].size());
	assert(m_offsets[sourceNumber][imageNumber] < m_originalPointProperties[(int)n][sourceNumber].size());
	return &(m_originalPointProperties[(int)n][sourceNumber][m_offsets[sourceNumber][imageNumber]]);
}

inline int ProjectedImagesInterface::getOriginalNumberOfTimeDelays(int sourceNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_originalTimeDelayInfo.size());
	return m_originalTimeDelayInfo[sourceNumber].size(); 
}

inline void ProjectedImagesInterface::getOriginalTimeDelay(int sourceNumber, int index, int *pImg, int *pPoint, float *pDelay) const	
{
	assert(pImg);
	assert(pPoint);
	assert(pDelay);
	assert(sourceNumber >= 0 && sourceNumber < m_originalTimeDelayInfo.size());
	assert(index >= 0 && index < m_originalTimeDelayInfo[sourceNumber].size());

	*pImg = m_originalTimeDelayInfo[sourceNumber][index].getImageIndex(); 
	*pPoint = m_originalTimeDelayInfo[sourceNumber][index].getPointIndex();
	*pDelay = m_originalTimeDelayInfo[sourceNumber][index].getTimeDelay(); 
}

} // end namespace

#endif // GRALE_PROJECTEDIMAGESINTERFACE_H

