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

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include "imagesdataextended.h"

#include "debugnew.h"

namespace grale
{

void ProjectedImagesInterface::storeOriginalData(const std::vector<ImagesDataExtended *> &images,
		                             bool storeOriginalIntensities, 
					     bool storeOriginalTimeDelays, 
					     bool storeOriginalShearInfo)
{
	m_numTotalPoints.resize(images.size());
	m_offsets.resize(images.size());
	m_numPoints.resize(images.size());
	m_distanceFractions.resize(images.size());

	for (int s = 0 ; s < images.size() ; s++)
	{
		ImagesDataExtended *pImgDat = images[s];
		int numImages = pImgDat->getNumberOfImages();
		int offset = 0;

		m_distanceFractions[s] = (float)(pImgDat->getDds()/pImgDat->getDs());

		m_offsets[s].resize(numImages);
		m_numPoints[s].resize(numImages);

		for (int i = 0 ; i < numImages ; i++)
		{
			int numPoints = pImgDat->getNumberOfImagePoints(i);

			m_offsets[s][i] = offset;
			m_numPoints[s][i] = numPoints;
			offset += numPoints;
		}

		m_numTotalPoints[s] = offset;
	}

	m_intensityScale = 1;

	m_originalIntensityFlags.resize(images.size());
	m_originalShearInfoFlags.resize(images.size());
	m_originalTimeDelayInfo.resize(images.size());

	for (int i = 0 ; i < images.size() ; i++)
	{
		m_originalShearInfoFlags[i] = false;
		m_originalIntensityFlags[i] = false;
	}

	if (storeOriginalIntensities)
	{
		double maxIntensity = 0;

		m_originalIntensities.resize(images.size());

		for (int source = 0 ; source < images.size() ; source++)
		{
			ImagesDataExtended *pImgDat = images[source];

			if (pImgDat->hasIntensities())
			{
				int numImages = pImgDat->getNumberOfImages();

				for (int i = 0 ; i < numImages ; i++)
				{
					int numPoints = pImgDat->getNumberOfImagePoints(i);

					for (int p = 0 ; p < numPoints ; p++)
					{
						double intens = ABS(pImgDat->getImagePointIntensity(i, p));

						if (intens > maxIntensity)
							maxIntensity = intens;
					}
				}
			}
		}

		if (maxIntensity == 0)
			m_intensityScale = 1;
		else
			m_intensityScale = maxIntensity/10.0;

		for (int source = 0 ; source < images.size() ; source++)
		{
			ImagesDataExtended *pImgDat = images[source];

			if (pImgDat->hasIntensities())
			{
				m_originalIntensityFlags[source] = true;
				m_originalIntensities[source].resize(m_numTotalPoints[source]);

				int numImages = pImgDat->getNumberOfImages();
				int point = 0;

				for (int i = 0 ; i < numImages ; i++)
				{
					int numPoints = pImgDat->getNumberOfImagePoints(i);

					for (int p = 0 ; p < numPoints ; p++, point++)
					{
						double intens = pImgDat->getImagePointIntensity(i, p)/m_intensityScale;

						m_originalIntensities[source][point] = (float)intens;
					}
				}
			}
		}
	}

	if (storeOriginalTimeDelays)
	{
		for (int source = 0 ; source < images.size() ; source++)
		{
			ImagesDataExtended *pImgDat = images[source];
			int numTimeDelays = pImgDat->getNumberOfTimeDelays();

			for (int i = 0 ; i < numTimeDelays ; i++)
			{
				int imgNumber, pointNumber;
				double timeDelay;

				pImgDat->getTimeDelay(i, &imgNumber, &pointNumber, &timeDelay);
				m_originalTimeDelayInfo[source].push_back(TimeDelayPoint(imgNumber, pointNumber, (float)timeDelay));
			}
		}	
	}

	if (storeOriginalShearInfo)
	{
		m_originalShearComponent1s.resize(images.size());
		m_originalShearComponent2s.resize(images.size());
		m_shearWeights.resize(images.size());

		for (int source = 0 ; source < images.size() ; source++)
		{
			ImagesDataExtended *pImgDat = images[source];

			if (pImgDat->hasShearInfo())
			{
				m_originalShearInfoFlags[source] = true;
				m_originalShearComponent2s[source].resize(m_numTotalPoints[source]);
				m_originalShearComponent1s[source].resize(m_numTotalPoints[source]);
				m_shearWeights[source].resize(m_numTotalPoints[source]);

				int numImages = pImgDat->getNumberOfImages();
				int point = 0;

				for (int i = 0 ; i < numImages ; i++)
				{
					int numPoints = pImgDat->getNumberOfImagePoints(i);

					for (int p = 0 ; p < numPoints ; p++, point++)
					{
						m_originalShearComponent1s[source][point] = (float)pImgDat->getShearComponent1(i, p);
						m_originalShearComponent2s[source][point] = (float)pImgDat->getShearComponent2(i, p);
						m_shearWeights[source][point] = (float)pImgDat->getShearWeight(i, p);
					}
				}
			}
		}
	}
}

/*
bool ProjectedImagesInterface::setDistanceFractions(const std::vector<float> &fractions)
{
	if (fractions.size() != m_distanceFractions.size())
		return false;

	m_distanceFractions = fractions;

	return true;
}
*/

} // end namespace
