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

#include "lensfitnesspointoverlap.h"
#include "fitnessutil.h"
#include <grale/vector2d.h>
#include <grale/polygon2d.h>
#include <grale/imagesdataextended.h>
#include <grale/triangle2d.h>
#include <grale/gridfunction.h>
#include <vector>
#include <iostream>
#include <limits>
#include <stdlib.h>

using namespace std;

namespace grale
{

LensFitnessPointOverlap::LensFitnessPointOverlap()
{
	m_initialized = false;
}

LensFitnessPointOverlap::~LensFitnessPointOverlap()
{
}

bool LensFitnessPointOverlap::init(double z_d, std::list<ImagesDataExtended *> &images,
		                           std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams)
{
	std::cerr << "NEW VERSION 2" << std::endl;

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

	m_deflectionFlags.resize(0);
	m_derivativeFlags.resize(0);
	m_potentialFlags.resize(0);
	m_distanceFractions.resize(0);
	m_sourceIndices.resize(0);
	
	std::list<ImagesDataExtended *>::const_iterator it;
	
	int imgIdx = 0;
	for (it = images.begin() ; it != images.end() ; it++, imgIdx++)
	{
		ImagesDataExtended *pImgDat = *it;
		int num;
		
		if ((num = pImgDat->getNumberOfImages()) <= 1)
		{
			setErrorString("An images data set doesn't contain multiple images");
			return false;
		}

		for (int i = 0 ; i < num ; i++)
		{
			if (pImgDat->getNumberOfImagePoints(i) != 1)
			{
				setErrorString("Each image must contain exactly one point");
				return false;
			}
		}

		m_deflectionFlags.push_back(true);
		m_derivativeFlags.push_back(false);
		m_potentialFlags.push_back(false);

		m_sourceIndices.push_back(imgIdx);

		double distFrac = pImgDat->getDds()/pImgDat->getDs();

		m_distanceFractions.push_back((float)distFrac);
	}

	m_initialized = true;
		
	return true;
}

bool LensFitnessPointOverlap::calculateOverallFitness(const ProjectedImagesInterface &interface, float *pFitnessValues) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}
	
	int numSources = interface.getNumberOfSources();
	float scale = getScaleFactor_PointImages(interface, m_sourceIndices, m_distanceFractions);
	float fitness = calculateOverlapFitness_PointImages(interface, m_sourceIndices, m_distanceFractions, scale);

	pFitnessValues[0] = fitness;
	return true;
}

string LensFitnessPointOverlap::getUsage() const
{
	return R"TEXT(
This module is kept only for backwards compatibility. 

Please use the 'general' module instead.

)TEXT";
}

LensFitnessPointOverlapNull::LensFitnessPointOverlapNull()
{
	m_initialized = false;
}

LensFitnessPointOverlapNull::~LensFitnessPointOverlapNull()
{
}

bool LensFitnessPointOverlapNull::init(double z_d, std::list<ImagesDataExtended *> &images, 
                                       std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams)
{
	std::cerr << "NEW VERSION 2" << std::endl;

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

	if (images.size()%2 != 0)
	{
		setErrorString("Number of images data sets should be a multiple of two");
		return false;
	}


	m_deflectionFlags.resize(0);
	m_derivativeFlags.resize(0);
	m_potentialFlags.resize(0);
	m_shortDeflectionFlags.resize(0);
	m_shortDerivativeFlags.resize(0);
	m_shortPotentialFlags.resize(0);
	m_distanceFractions.resize(0);
	m_nullWeights.resize(0);
	m_sourceIndices.resize(0);
	m_nullIndices.resize(0);
	m_shortSourceIndices.resize(0);
	
	int numSources = images.size()/2;
	m_nullTriangles.resize(numSources);

	std::list<ImagesDataExtended *>::const_iterator it;
	int count = 0;
	int multipleImageSources = 0;
	
	for (it = images.begin() ; it != images.end() ; it++, count++)
	{
		ImagesDataExtended *pImgDat = *it;
		int sourcePos = count/2;

		if (( count % 2 ) == 0)
		{
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

			m_deflectionFlags.push_back(true);
			m_derivativeFlags.push_back(false);
			m_potentialFlags.push_back(false);
			m_shortDeflectionFlags.push_back(true);
			m_shortDerivativeFlags.push_back(false);
			m_shortPotentialFlags.push_back(false);

			m_shortSourceIndices.push_back(sourcePos);
			m_sourceIndices.push_back(count);

			shortImages.push_back(pImgDat);

			double distFrac = pImgDat->getDds()/pImgDat->getDs();

			m_distanceFractions.push_back((float)distFrac);
	
			string errStr;
			shared_ptr<ImagesDataExtended> grpImg = addGroupsToPointImages(*pImgDat, errStr);
			if (!grpImg.get())
			{
				setErrorString(errStr);
				return false;
			}
			m_pointGroups.add(grpImg.get());
		}
		else // null space data
		{
			float nullWeight = 1.0f; // default weight is 1, reproduces old behaviour

			m_deflectionFlags.push_back(true);
			m_derivativeFlags.push_back(false);
			m_potentialFlags.push_back(false);

			m_nullIndices.push_back(count);

			int num;
			
			if ((num = pImgDat->getNumberOfImages()) != 1)
			{
				setErrorString("Null space data should contain only one image");
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
					setErrorString("Unable to get extra parameter specifying null space weight: " + pImgDat->getErrorString());
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

				std::vector<TriangleIndices>::const_iterator triangIt;
				int i;

				for (i = 0, triangIt = t.begin() ; triangIt != t.end() ; triangIt++, i++)
					m_nullTriangles[sourcePos][i] = *triangIt;
			}
		}
	}

	if (multipleImageSources < 1)
	{
		setErrorString("No source with multiple images is present");
		return false;
	}

	m_initialized = true;
		
	return true;
}

bool LensFitnessPointOverlapNull::calculateMassScaleFitness(const ProjectedImagesInterface &interface, float &fitness) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return true;
	}
	
	int numSources = interface.getNumberOfSources();
	float scale = getScaleFactor_PointImages(interface, m_shortSourceIndices, m_distanceFractions);
	fitness = calculateOverlapFitness_PointImages(interface, m_shortSourceIndices, m_distanceFractions, scale);

	return true;
}

bool LensFitnessPointOverlapNull::calculateOverallFitness(const ProjectedImagesInterface &interface, float *pFitnessValues) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}

	int numSources = interface.getNumberOfSources();
	float scale = getScaleFactor_PointImages(interface, m_sourceIndices, m_distanceFractions);

	pFitnessValues[0] = calculateOverlapFitness_PointImages(interface, m_sourceIndices, m_distanceFractions, scale);
	pFitnessValues[1] = calculateNullFitness_PointImages(m_pointGroups, interface, m_sourceIndices, m_nullIndices, m_nullTriangles, m_nullWeights);

	return true;
}

string LensFitnessPointOverlapNull::getUsage() const
{
	return R"TEXT(
This module is kept only for backwards compatibility. 

Please use the 'general' module instead.

)TEXT";
}

} // end namespace

