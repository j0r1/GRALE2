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

namespace grale
{

void ProjectedImagesInterface::storeOriginalData(const std::vector<ImagesDataExtended *> &images)
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

	m_originalTimeDelayInfo.resize(images.size());

	m_originalPropertyFlags.resize(ImagesData::MaxProperty);
	for (auto &p : m_originalPropertyFlags)
		p.resize(images.size(), false);

	m_originalPointProperties.resize(ImagesData::MaxProperty);
	for (size_t prop = 0 ; prop < m_originalPointProperties.size() ; prop++)
	{
		ImagesData::PropertyName propName = (ImagesData::PropertyName)prop;
		auto &pp = m_originalPointProperties[prop];
		pp.resize(images.size());

		for (size_t source = 0 ; source < images.size() ; source++)
		{
			ImagesDataExtended *pImgDat = images[source];
			if (pImgDat->hasProperty(propName))
			{
				m_originalPropertyFlags[prop][source] = true;
				pp[source].resize(m_numTotalPoints[source]);

				auto &ppSrc = pp[source];
				int numImages = pImgDat->getNumberOfImages();
				int point = 0;

				for (int i = 0 ; i < numImages ; i++)
				{
					int numPoints = pImgDat->getNumberOfImagePoints(i);
					for (int p = 0 ; p < numPoints ; p++, point++)
						ppSrc[point] = (float)pImgDat->getImagePointProperty(propName, i, p);
				}
			}
		}
	}

	if (1)
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
}

} // end namespace
