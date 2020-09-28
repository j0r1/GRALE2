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
	int getNumberOfSources() const;

	/** Returns the number of images stored for a specific source. */
	int getNumberOfImages(int sourceNumber) const;
	
	/** Returns the total number of points of all images corresponding to a source index. */
	int getNumberOfImagePoints(int sourceNumber) const;

	/** Returns the number of points stored for a particular image of a particular source. */
	int getNumberOfImagePoints(int sourceNumber, int imageNumber) const;

	/** Returns true if measured intensity information was stored for a specific source index. */
	bool hasOriginalIntensities(int sourceNumber) const;

	/** Returns true if measured shear info was stored for the points of a specific source. */
	bool hasOriginalShearInfo(int sourceNumber) const;

	/** Returns the measured intensities stored for all images of a specific source number
	 *  (expressed in terms of the intensity scale, see getIntensityScale). */
	const float *getOriginalIntensities(int sourceNumber) const;

	/** Returns the measured intensities stored for a specific image of a specific source.
	 *  (expressed in terms of the intensity scale, see getIntensityScale). */
	const float *getOriginalIntensities(int sourceNumber, int imageNumber) const;

	/** Returns the stored gamma1 values for all images of a specific source. */
	const float *getOriginalShearComponent1s(int sourceNumber) const;

	/** Returns the stored gamma1 values for a specific image of a specific source. */
	const float *getOriginalShearComponent1s(int sourceNumber, int imageNumber) const;

	/** Returns the stored gamma2 values for all images of a specific source. */
	const float *getOriginalShearComponent2s(int sourceNumber) const;

	/** Returns the stored gamma2 values for a specific image of a specific source. */
	const float *getOriginalShearComponent2s(int sourceNumber, int imageNumber) const;

	const float *getOriginalShearWeights(int sourceNumber) const;
	const float *getOriginalShearWeights(int sourceNumber, int imageNumber) const;

	/** Returns the number of time delays that were stored for a specific source index. */
	int getOriginalNumberOfTimeDelays(int sourceNumber) const;

	/** For a specific source and specific index (from 0 to getOriginalNumberOfTimeDelays - 1),
	 *  the image index for a time delay value will be stored in \c pImg, the point index in \c pPoint
	 *  and the actual measured time delay in \c pDelay. */
	void getOriginalTimeDelay(int sourceNumber, int index, int *pImg, int *pPoint, float *pDelay) const;

	/** Returns the scale in which stored and calculated image point intensities are expressed. */
	double getIntensityScale() const;

	// For internal use only
	//bool setDistanceFractions(const std::vector<float> &fractions);

	float getDistanceFraction(int sourcenum) const { return m_distanceFractions[sourcenum]; }

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

	void dump(bool magnifyFlux) const
	{
		int numSources = getNumberOfSources();

		for (int s = 0 ; s < numSources ; s++)
		{
			int numImages = getNumberOfImages(s);

			for (int i = 0 ; i < numImages ; i++)
			{
				int numPoints = getNumberOfImagePoints(s, i);

				for (int p = 0 ; p < numPoints ; p++)
				{
					Vector2D<float> beta = getBetas(s, i)[p];
					Vector2D<double> beta2(beta.getX(), beta.getY());
					double intens = getOriginalIntensities(s, i)[p];

					if (magnifyFlux)
						intens *= getInverseMagnifications(s, i)[p];

					beta2 *= getAngularScale();
					intens *= getIntensityScale();

					printf("%10.10g %10.10g %10.10g\n",(double)(beta2.getX()/ANGLE_ARCSEC),(double)(beta2.getY()/ANGLE_ARCSEC),intens);
				}
			}
		}
	}
protected:
	void storeOriginalData(const std::vector<ImagesDataExtended *> &images,
			       bool storeOriginalIntensities, bool storeOriginalTimeDelays, bool storeOriginalShearInfo);

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

	double m_intensityScale;
	std::vector<float> m_distanceFractions;
	std::vector<bool> m_originalIntensityFlags;
	std::vector<bool> m_originalShearInfoFlags;
	std::vector<std::vector<float> > m_originalIntensities;
	std::vector<std::vector<float> > m_originalShearComponent1s;
	std::vector<std::vector<float> > m_originalShearComponent2s;
	std::vector<std::vector<float> > m_shearWeights;
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

inline bool ProjectedImagesInterface::hasOriginalIntensities(int sourceNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_originalIntensityFlags.size());
	return m_originalIntensityFlags[sourceNumber]; 
}

inline bool ProjectedImagesInterface::hasOriginalShearInfo(int sourceNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_originalShearInfoFlags.size());
	return m_originalShearInfoFlags[sourceNumber]; 
}

inline const float *ProjectedImagesInterface::getOriginalIntensities(int sourceNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_originalIntensities.size());
	assert(m_originalIntensities[sourceNumber].size() > 0);
	return &(m_originalIntensities[sourceNumber][0]); 
}

inline const float *ProjectedImagesInterface::getOriginalIntensities(int sourceNumber, int imageNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_originalIntensities.size());
	assert(sourceNumber < m_offsets.size());
	assert(imageNumber >= 0 && imageNumber < m_offsets[sourceNumber].size());
	assert(m_offsets[sourceNumber][imageNumber] < m_originalIntensities[sourceNumber].size());
	return &(m_originalIntensities[sourceNumber][m_offsets[sourceNumber][imageNumber]]); 
}

inline const float *ProjectedImagesInterface::getOriginalShearComponent1s(int sourceNumber) const
{ 
	assert(sourceNumber >= 0 && sourceNumber < m_originalShearComponent1s.size());
	assert(m_originalShearComponent1s[sourceNumber].size() > 0);
	return &(m_originalShearComponent1s[sourceNumber][0]); 
}

inline const float *ProjectedImagesInterface::getOriginalShearComponent1s(int sourceNumber, int imageNumber) const
{ 
	assert(sourceNumber >= 0 && sourceNumber < m_originalShearComponent1s.size());
	assert(sourceNumber < m_offsets.size());
	assert(imageNumber >= 0 && imageNumber < m_offsets[sourceNumber].size());
	assert(m_offsets[sourceNumber][imageNumber] < m_originalShearComponent1s[sourceNumber].size());
	return &(m_originalShearComponent1s[sourceNumber][m_offsets[sourceNumber][imageNumber]]); 
}

inline const float *ProjectedImagesInterface::getOriginalShearComponent2s(int sourceNumber) const
{ 
	assert(sourceNumber >= 0 && sourceNumber < m_originalShearComponent2s.size());
	assert(m_originalShearComponent2s[sourceNumber].size() > 0);
	return &(m_originalShearComponent2s[sourceNumber][0]); 
}

inline const float *ProjectedImagesInterface::getOriginalShearComponent2s(int sourceNumber, int imageNumber) const
{ 
	assert(sourceNumber >= 0 && sourceNumber < m_originalShearComponent2s.size());
	assert(sourceNumber < m_offsets.size());
	assert(imageNumber >= 0 && imageNumber < m_offsets[sourceNumber].size());
	assert(m_offsets[sourceNumber][imageNumber] < m_originalShearComponent2s[sourceNumber].size());
	return &(m_originalShearComponent2s[sourceNumber][m_offsets[sourceNumber][imageNumber]]); 
}

inline const float *ProjectedImagesInterface::getOriginalShearWeights(int sourceNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_shearWeights.size());
	assert(m_shearWeights[sourceNumber].size() > 0);
	return &(m_shearWeights[sourceNumber][0]);
}

inline const float *ProjectedImagesInterface::getOriginalShearWeights(int sourceNumber, int imageNumber) const
{
	assert(sourceNumber >= 0 && sourceNumber < m_shearWeights.size());
	assert(sourceNumber < m_offsets.size());
	assert(imageNumber >= 0 && imageNumber < m_offsets[sourceNumber].size());
	assert(m_offsets[sourceNumber][imageNumber] < m_shearWeights[sourceNumber].size());
	return &(m_shearWeights[sourceNumber][m_offsets[sourceNumber][imageNumber]]);
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

inline double ProjectedImagesInterface::getIntensityScale() const
{ 
	return m_intensityScale;
}

} // end namespace

#endif // GRALE_PROJECTEDIMAGESINTERFACE_H

