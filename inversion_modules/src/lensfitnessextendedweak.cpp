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

#include "lensfitnessextendedweak.h"
#include <grale/vector2d.h>
#include <grale/polygon2d.h>
#include <grale/imagesdataextended.h>
#include <grale/triangle2d.h>
#include <grale/gridfunction.h>
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;

namespace grale
{

LensFitnessExtendedWeak::LensFitnessExtendedWeak()
{
	m_initialized = false;
}

LensFitnessExtendedWeak::~LensFitnessExtendedWeak()
{
}

bool LensFitnessExtendedWeak::init(double z_d, std::list<ImagesDataExtended *> &images, 
                                   std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams)
{
	if (m_initialized)
	{
		setErrorString("Already initialized");
		return false;
	}

	m_pointGroups.clear();

	std::list<ImagesDataExtended *>::const_iterator it, weakStartPos;
	int count = 0;
	
	if (images.size() == 0)
	{
		setErrorString("No available images data");
		return false;
	}

	m_weakStartPos = images.size(); // init at the end of the images, in case there isn't any weak lensing info
	weakStartPos = images.end();

	// Find the start of the weak lensing info
	for (it = images.begin(), count = 0 ; it != images.end() ; it++, count++)
	{
		ImagesDataExtended *pImgDat = *it;

		if (count%2 == 0)
		{
			if (pImgDat->getNumberOfExtraParameters() == 3)
			{
				double val = 0;
				if (pImgDat->getExtraParameter("0", val))
				{
					if (val < 0) // use something negative to mark the start of weak lensing data
					{
						m_weakStartPos = count;
						weakStartPos = it;
						break; // ok, we found it
					}
				}
			}
		}
	}

	if (m_weakStartPos%2 != 0)
	{
		setErrorString("Expected an even number of image data sets for the strong lensing data, each time extended images and null space");
		return false;
	}

	int numSources = m_weakStartPos/2; // everything before the start of the weak lensing data is a 'normal' source

	m_nullTriangles.resize(numSources);
	m_nullTriangleAreas.resize(numSources);
	m_nullWeights.resize(numSources);
	
	m_sourceIndices.clear();
	m_shortSourceIndices.clear();
	m_nullIndices.clear();
	m_useRectFlags.clear();
	m_usePointGroupFlags.clear();
	m_weakIndices.clear();

	m_shortNumberOfSources = numSources;
	
	for (it = images.begin(), count = 0 ; it != images.end() && count < m_weakStartPos ; it++, count++)
	{
		int sourcePos = count/2;
		ImagesDataExtended *pImgDat = *it;
		
		m_totalDeflectionFlags.push_back(true);
		m_weakOneFlags.push_back(false);

		if (count%2 == 0) // actual images
		{
			m_pointGroups.add(pImgDat);

			shortImages.push_back(pImgDat);
			m_shortSourceIndices.push_back(sourcePos);
			m_sourceIndices.push_back(count);
			m_useRectFlags.push_back(true);
			m_usePointGroupFlags.push_back(true);

			m_shortDeflectionFlags.push_back(true);
			m_shortDerivativeFlags.push_back(false);
			m_shortPotentialFlags.push_back(false);

			m_totalDerivativeFlags.push_back(false);
			m_totalPotentialFlags.push_back(false);

			int num;
			
			if ((num = pImgDat->getNumberOfImages()) < 1)
			{
				setErrorString("An images data set doesn't contain any images");
				return false;
			}

			if (num == 1)
			{
				if (pImgDat->getNumberOfImagePoints(0) != 1)
				{
					setErrorString("In case only a single image is present, it must be a single point (used for null space)");
					return false;
				}
			}
			else
			{
				for (int i = 0 ; i < num ; i++)
				{
					if (pImgDat->getNumberOfImagePoints(i) <= 2)
					{
						setErrorString("Each image must contain at least three points");
						return false;
					}
				}
			}

			int numExtraParams = pImgDat->getNumberOfExtraParameters();

			if (numExtraParams != 0)
			{
				setErrorString("No extra parameters are supported for image data");
				return false;
			}
		}
		else // null space data
		{
			float nullWeight = 1.0f; // default weight is 1, reproduces old behaviour

			m_nullIndices.push_back(count);

			m_totalDerivativeFlags.push_back(false);
			m_totalPotentialFlags.push_back(false);

			if (pImgDat->getNumberOfImages() != 1)
			{
				setErrorString("Null space data should consist of only one image");
				return false;
			}

			const int numExtra = pImgDat->getNumberOfExtraParameters();
			if (numExtra > 1)
			{
				setErrorString("Only one extra parameter specifying a weigth is allowed for null space data");
				return false;
			}

			if (numExtra == 1)
			{
				double val = -1;

				if (!pImgDat->getExtraParameter("0", val))
				{
					setErrorString("Can't get null space weight: " + pImgDat->getErrorString());
					return false;
				}

				if (val < 0 || val >= (double)std::numeric_limits<float>::max())
				{
					setErrorString("Extra parameter specifying weight of null space is either negative or too large");
					return false;
				}

				nullWeight = (float)val;
			}

			m_nullWeights[sourcePos] = nullWeight;

			if (pImgDat->getNumberOfImagePoints(0) == 1) // we use this to signal that the null space should be skipped
			{
				m_nullTriangles[sourcePos].resize(0);
				m_nullTriangleAreas[sourcePos].resize(0);
			}
			else // real null space
			{
				if (!pImgDat->hasTriangulation())
				{
					setErrorString("Null space data doesn't contain a triangulation");
					return false;
				}

				std::vector<TriangleIndices> t;

				if (!pImgDat->getTriangles(0, t))
				{
					setErrorString(std::string("Unable to obtain triangulation data from the null space data: ") + (*it)->getErrorString());
					return false;
				}
				
				m_nullTriangles[sourcePos].resize(t.size());
				m_nullTriangleAreas[sourcePos].resize(t.size());

				std::vector<TriangleIndices>::const_iterator triangIt;
				std::vector<Vector2D<double> > points;
				int i;

				for (i = 0, triangIt = t.begin() ; triangIt != t.end() ; triangIt++, i++)
					m_nullTriangles[sourcePos][i] = (*triangIt);
				
				points.resize(pImgDat->getNumberOfImagePoints(0));

				for (i = 0 ; i < points.size() ; i++)
					points[i] = pImgDat->getImagePointPosition(0, i);
				
				for (i = 0 ; i < m_nullTriangles[sourcePos].size() ; i++)
				{
					Vector2D<double> p1 = points[m_nullTriangles[sourcePos][i].getIndex(0)];
					Vector2D<double> p2 = points[m_nullTriangles[sourcePos][i].getIndex(1)];
					Vector2D<double> p3 = points[m_nullTriangles[sourcePos][i].getIndex(2)];

					Triangle2D<double> t(p1,p2,p3);
					m_nullTriangleAreas[sourcePos][i] = t.getArea()/(ANGLE_ARCSEC*ANGLE_ARCSEC);
				}
			}
		}
	}

	for (it = weakStartPos ; it != images.end() ; it++, count++)
	{
		ImagesDataExtended *pImgDat = *it;

		m_weakIndices.push_back(count);
		m_totalDeflectionFlags.push_back(false);
		m_totalDerivativeFlags.push_back(true);
		m_totalPotentialFlags.push_back(false);
		m_weakOneFlags.push_back(true);

		if (pImgDat->getNumberOfImages() != 1)
		{
			setErrorString("Each images data instance can only contain one 'image'");
			return false;
		}

		if (!pImgDat->hasShearInfo())
		{
			setErrorString("No shear info is present in a data set");
			return false;
		}

		if (pImgDat->getNumberOfExtraParameters() != 3)
		{
			setErrorString("Weak lensing points need three extra parameters: any negative value (as a marker), value 1 or 2 for reduced shear or regular shear, and a threshold for |1-kappa|");
			return false;
		}

		double marker, shearType, threshold;

		if (!pImgDat->getExtraParameter("0", marker) || !pImgDat->getExtraParameter("1", shearType) ||
		    !pImgDat->getExtraParameter("2", threshold))
		{
			setErrorString("Unable to get marker, shear type or threshold parameter: " + pImgDat->getErrorString());
			return false;
		}

		if (marker >= 0)
		{
			setErrorString("Weak lensing data must have a negative value as a first extra parameter");
			return false;
		}

		if (shearType == 1.0)
			m_reducedShear.push_back(RealReducedShear);
		else if (shearType == 2.0)
			m_reducedShear.push_back(RealShear);
		else
		{
			setErrorString("As a second extra parameter, weak lensing data must have either 1 (reduced shear) or 2 (regular shear)");
			return false;
		}

		if (threshold < 0)
		{
			setErrorString("The threshold value for |1-kappa| must be positive or zero");
			return false;
		}

		m_thresholds.push_back(threshold);
	}
		
	std::cerr << "m_weakStartPos = " << m_weakStartPos << std::endl;
	std::cerr << "Strong images (with null space) = " << numSources << std::endl;
	std::cerr << "Weak lensing image groups = " << (images.size()-m_weakStartPos) << std::endl;
	
	m_initialized = true;

	return true;
}

bool LensFitnessExtendedWeak::calculateMassScaleFitness(const ProjectedImagesInterface &iface, float &fitness) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}
	
	int numsources = iface.getNumberOfSources();
	if (numsources != m_shortNumberOfSources)
	{
		setErrorString("Wrong number of sources!");
		return false;
	}

	fitness = calculateOverlapFitness_Extended(m_pointGroups, iface, m_shortSourceIndices, m_useRectFlags, m_usePointGroupFlags);
	return true;
}

bool LensFitnessExtendedWeak::calculateOverallFitness(const ProjectedImagesInterface &iface, float *fitnessvalues) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}

	float posfitness = calculateOverlapFitness_Extended(m_pointGroups, iface, m_sourceIndices, m_useRectFlags,
			                                            m_usePointGroupFlags);
	float nullfitness = calculateNullFitness_ExtendedImages(iface, m_sourceIndices, m_nullIndices, m_nullTriangles, 
			                                                m_nullTriangleAreas, m_nullWeights);

	fitnessvalues[0] = posfitness;
	fitnessvalues[1] = nullfitness;

	// Now we do the weak lensing part

	int endPos = iface.getNumberOfSources();
	fitnessvalues[2] = calculateWeakLensingFitness(iface, m_weakIndices, m_reducedShear, m_thresholds);

	return true;
}

string LensFitnessExtendedWeak::getUsage() const
{
	return R"TEXT(
This module is kept only for backwards compatibility. 

Please use the 'general' module instead.

)TEXT";
}

} // end namespace

