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

#include "lensfitnessweak.h"
#include "fitnessutil.h"
#include <grale/imagesdataextended.h>

using namespace std;

namespace grale
{

LensFitnessWeak::LensFitnessWeak()
{
	m_initialized = false;
}

LensFitnessWeak::~LensFitnessWeak()
{
}

bool LensFitnessWeak::init(double z_d, std::list<ImagesDataExtended *> &images, 
                           std::list<ImagesDataExtended *> &massScaleImages, const ConfigurationParameters *pParams)
{
	if (m_initialized)
	{
		setErrorString("Already initialized");
		return false;
	}

	std::list<ImagesDataExtended *>::const_iterator it;

	m_oneVector.resize(0);
	m_zeroVector.resize(0);
	m_weakIndices.resize(0);

	int count = 0;
	for (it = images.begin() ; it != images.end() ; it++, count++)
	{
		ImagesDataExtended *pImgDat = *it;

		if (pImgDat->getNumberOfImages() != 1)
		{
			setErrorString("Each images data instance can only contain one image");
			return false;
		}

		if (!pImgDat->hasShearInfo())
		{
			setErrorString("No shear info is present");
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

		m_weakIndices.push_back(count);
		m_thresholds.push_back(threshold);

		m_oneVector.push_back(1);
		m_zeroVector.push_back(0);
	}

	m_initialized = true;
	return true;
}

bool LensFitnessWeak::calculateOverallFitness(const ProjectedImagesInterface &interface, float *pFitnessValues) const
{
	if (!m_initialized)
	{
		setErrorString("Not initialized");
		return false;
	}

	int numSources = interface.getNumberOfSources();

	pFitnessValues[0] = calculateWeakLensingFitness(interface, m_weakIndices, m_reducedShear, m_thresholds);
	return true;
}

string LensFitnessWeak::getUsage() const
{
	return R"TEXT(
This module is kept only for backwards compatibility. 

Please use the 'general' module instead.

)TEXT";
}

} // end namespace
