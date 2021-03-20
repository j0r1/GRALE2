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

#include "lensfitnesssimplerectangles.h"
#include <grale/fitnessutil.h>
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

bool LensFitnessSimpleRectangles::init(double z_d, std::list<ImagesDataExtended *> &images, 
                                       std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams)
{
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
	
	std::list<ImagesDataExtended *>::const_iterator it;
	
	int count = 0;
	for (it = images.begin() ; it != images.end() ; it++, count++)
	{
		int num;
		
		if ((num = (*it)->getNumberOfImages()) <= 1)
		{
			setErrorString("An images data set doesn't contain multiple images");
			return false;
		}

		for (int i = 0 ; i < num ; i++)
		{
			if ((*it)->getNumberOfImagePoints(i) <= 2)
			{
				setErrorString("Each image must contain at least three points");
				return false;
			}
		}

		m_sourceIndices.push_back(count);
		m_useRectFlags.push_back(true);
		m_usePointGroupFlags.push_back(true);

		m_deflectionFlags.push_back(true);
		m_derivativeFlags.push_back(false);
		m_potentialFlags.push_back(false);
	}

	m_pointGroupInfo.init(images);
	m_initialized = true;
		
	return true;
}

bool LensFitnessSimpleRectangles::calculateMassScaleFitness(const ProjectedImagesInterface &inf, float &fitness) const
{
	int numsources = inf.getNumberOfSources();
	fitness = calculateOverlapFitness_Extended(m_pointGroupInfo, inf, m_sourceIndices, m_useRectFlags, m_usePointGroupFlags);
	return true;
}

bool LensFitnessSimpleRectangles::calculateOverallFitness(const ProjectedImagesInterface &iface, float *fitnessvalues) const
{
	if (!m_initialized)
	{
		setErrorString("Not initialized");
		return false;
	}

	int numsources = iface.getNumberOfSources();

	fitnessvalues[0] = calculateOverlapFitness_Extended(m_pointGroupInfo, iface, m_sourceIndices, m_useRectFlags, 
			                                            m_usePointGroupFlags);

	return true;
}

string LensFitnessSimpleRectangles::getUsage() const
{
	return R"TEXT(
This module is kept only for backwards compatibility. 

Please use the 'general' module instead.

)TEXT";
}

} // end namespace

