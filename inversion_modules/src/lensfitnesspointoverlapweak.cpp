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

#include "lensfitnesspointoverlapweak.h"
#include "fitnessutil.h"
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

LensFitnessPointOverlapNullWeak::LensFitnessPointOverlapNullWeak()
{
	m_initialized = false;
}

LensFitnessPointOverlapNullWeak::~LensFitnessPointOverlapNullWeak()
{
}

bool LensFitnessPointOverlapNullWeak::init(double z_d, std::list<ImagesDataExtended *> &images, 
                                           std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams)
{
	std::cerr << "STRONG & WEAK TEST" << std::endl;

	if (m_initialized)
	{
		setErrorString("Already initialized");
		return false;
	}
	
	if (images.size() == 0)
	{
		setErrorString("No available images data");
		return false;
	}

	//if (images.size()%2 != 0)
	//{
	//	setErrorString("Number of images data sets should be a multiple of two");
	//	return false;
	//}

	m_deflectionFlags.resize(0);
	m_derivativeFlags.resize(0);
	m_potentialFlags.resize(0);
	m_shortDeflectionFlags.resize(0);
	m_shortDerivativeFlags.resize(0);
	m_shortPotentialFlags.resize(0);
	m_distanceFractions.resize(0);
	m_weakOneFlags.resize(0);
	m_nullWeights.resize(0);

	m_shortSourceIndices.clear();
	m_sourceIndices.clear();
	m_nullIndices.clear();
	m_weakIndices.clear();
	
	// This will allocate too much room, but that doesn't really matter
	int numSources = images.size()/2;
	m_nullTriangles.resize(numSources);

	std::list<ImagesDataExtended *>::const_iterator it, weakStartPos;
	int count = 0;
	int multipleImageSources = 0;
	
	weakStartPos = images.end(); // just an initialization
	m_weakStartPos = images.size();

	for (it = images.begin() ; it != images.end() ; it++, count++)
	{
		ImagesDataExtended *pImgDat = *it;
		int sourcePos = count/2;

		if (( count % 2 ) == 0)
		{
			// Check if we're starting with weak lensing data
			if (pImgDat->getNumberOfExtraParameters() == 3)
			{
				double val = 0;
				if (pImgDat->getExtraParameter("0", val))
				{
					if (val < 0) // use something negative to mark the start of weak lensing data
					{
						m_weakStartPos = count;
						weakStartPos = it;
						break; // quit this loop and continue with the weak lensing handling
					}
				}
			}

			// Just strong lensing, continue 
			int num;
			
			if ((num = pImgDat->getNumberOfImages()) < 1)
			{
				setErrorString("An images data set doesn't contain any images");
				return false;
			}

			if (num > 1)
				multipleImageSources++;

			for (int i = 0 ; i < num ; i++)
			{
				if (pImgDat->getNumberOfImagePoints(i) != 1)
				{
					setErrorString("Each image must contain exactly one point");
					return false;
				}
			}

			m_sourceIndices.push_back(count);
			m_shortSourceIndices.push_back(sourcePos);

			m_deflectionFlags.push_back(true);
			m_derivativeFlags.push_back(false);
			m_potentialFlags.push_back(false);
			m_shortDeflectionFlags.push_back(true);
			m_shortDerivativeFlags.push_back(false);
			m_shortPotentialFlags.push_back(false);
			m_weakOneFlags.push_back(false);

			shortImages.push_back(pImgDat);

			double distFrac = pImgDat->getDds()/pImgDat->getDs();

			m_distanceFractions.push_back((float)distFrac);
		}
		else // null space data
		{
			float nullWeight = 1.0f; // default weight is 1, reproduces old behaviour

			m_deflectionFlags.push_back(true);
			m_derivativeFlags.push_back(false);
			m_potentialFlags.push_back(false);
			m_weakOneFlags.push_back(false);

			m_nullIndices.push_back(count);

			int num;
			
			if ((num = pImgDat->getNumberOfImages()) != 1)
			{
				setErrorString("Null space data should contain only one image");
				return false;
			}
			if (!pImgDat->hasTriangulation())
			{
				setErrorString("Null space data doesn't contain a triangulation");
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

			m_nullWeights.push_back(nullWeight);

			if (pImgDat->getNumberOfImagePoints(0) == 1) // we use this to signal that the null space should be skipped
			{
				m_nullTriangles[sourcePos].resize(0);
			}
			else // real null space
			{
				std::vector<TriangleIndices> t;

				if (!pImgDat->getTriangles(0, t))
				{
					setErrorString(std::string("Unable to obtain triangulation data from the null space data: ") + (*it)->getErrorString());
					return false;
				}
				
				m_nullTriangles[sourcePos].resize(t.size());

				std::vector<TriangleIndices>::const_iterator triangIt;
				int i;

				for (i = 0, triangIt = t.begin() ; triangIt != t.end() ; triangIt++, i++)
					m_nullTriangles[sourcePos][i] = *triangIt;
			}
		}
	}

	for (it = weakStartPos ; it != images.end() ; it++, count++)
	{
		ImagesDataExtended *pImgDat = *it;

		m_deflectionFlags.push_back(false);
		m_derivativeFlags.push_back(true);
		m_potentialFlags.push_back(false);
		m_weakOneFlags.push_back(true);
		m_weakIndices.push_back(count);

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
			m_reducedShear.push_back(true);
		else if (shearType == 2.0)
			m_reducedShear.push_back(false);
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

	if (multipleImageSources < 1)
	{
		setErrorString("No source with multiple images is present");
		return false;
	}

	m_initialized = true;
		
	std::cerr << "m_weakStartPos = " << m_weakStartPos << std::endl;
	std::cerr << "Strong images (with null space) = " << m_weakStartPos/2 << std::endl;
	std::cerr << "Weak lensing image groups = " << (images.size()-m_weakStartPos) << std::endl;

	return true;
}

bool LensFitnessPointOverlapNullWeak::calculateMassScaleFitness(const ProjectedImagesInterface &interface, float &fitness) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}
	
	// for the 'short' version, the interface just contains the strong lensing images
	int numSources = interface.getNumberOfSources();

	float scale = getScaleFactor_PointImages(interface, m_shortSourceIndices, m_distanceFractions);
	fitness = calculateOverlapFitness_PointImages(interface, m_shortSourceIndices, m_distanceFractions, scale);

	return true;
}

bool LensFitnessPointOverlapNullWeak::calculateOverallFitness(const ProjectedImagesInterface &interface, float *pFitnessValues) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}

	//int numSources = interface.getNumberOfSources();
	int numSources = m_weakStartPos; // start of weak lensing data marks the end of strong lensing data. numSources here is just used for the strong lensing stuff
	float nullFitness = 0;

	float scale = getScaleFactor_PointImages(interface, m_sourceIndices, m_distanceFractions);

	pFitnessValues[0] = calculateOverlapFitness_PointImages(interface, m_sourceIndices, m_distanceFractions, scale);
	pFitnessValues[1] = calculateNullFitness_PointImages(interface, m_sourceIndices, m_nullIndices, m_nullTriangles, m_nullWeights);

	// Now we do the weak lensing part
	int endPos = interface.getNumberOfSources();
	pFitnessValues[2] = calculateWeakLensingFitness(interface, m_weakIndices, m_reducedShear, m_thresholds);

	return true;
}

string LensFitnessPointOverlapNullWeak::getUsage() const
{
	return R"TEXT(
This module is kept only for backwards compatibility. 

Please use the 'general' module instead.

)TEXT";
}

} // end namespace

