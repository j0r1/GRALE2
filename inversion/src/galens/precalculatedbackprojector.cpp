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
#include "precalculatedbackprojector.h"
#include "imagesdataextended.h"
#include "constants.h"
#include <memory>

namespace grale
{

using namespace std;

PreCalculatedBackProjector::PreCalculatedBackProjector()
{
}

PreCalculatedBackProjector::~PreCalculatedBackProjector()
{
}

bool PreCalculatedBackProjector::init(const std::vector<ImagesData *> &images, 
		                              const std::vector<ImagesData *> &correspondingSources)
{
	// Keep namespace clean
	{
		vector<ImagesDataExtended *> imagesVector;
		vector<unique_ptr<ImagesDataExtended>> imagesVectorAutoDelete;

		for (auto it : images)
			imagesVectorAutoDelete.push_back(make_unique<ImagesDataExtended>(*it));

		for (auto &it : imagesVectorAutoDelete)
			imagesVector.push_back(it.get());

		storeOriginalData(imagesVector);
	}

	if (images.size() == 0)
	{
		setErrorString("No images data sets were specified");
		return false;
	}

	if (images.size() != correspondingSources.size())
	{
		setErrorString("There should be an equal amount of data sets for images and sources");
		return false;
	}

	m_betas.resize(images.size());
	m_thetas.resize(images.size());

	double maxx = -numeric_limits<double>::max(), minx = numeric_limits<double>::max();
    double maxy = -numeric_limits<double>::max(), miny = numeric_limits<double>::max();
	
	// Determine the angular scale and store some stuff

	for (auto pImgDat : images)
	{
		assert(pImgDat);
		int numImages = pImgDat->getNumberOfImages();
		
		for (int j = 0 ; j < numImages ; j++)
		{
			int numPoints = pImgDat->getNumberOfImagePoints(j);

			for (int k = 0 ; k < numPoints ; k++)
			{
				Vector2Dd theta = pImgDat->getImagePointPosition(j, k);
				
				minx = MIN(minx,theta.getX());
				maxx = MAX(maxx,theta.getX());
				miny = MIN(miny,theta.getY());
				maxy = MAX(maxy,theta.getY());
			}			
		}
	}
	
	double dx = maxx-minx;
	double dy = maxy-miny;
	m_angularScale = sqrt(dx*dx+dy*dy);
	
	for (size_t s = 0 ; s < images.size() ; s++)
	{
		ImagesData *pImg = images[s];
		ImagesData *pSrc = correspondingSources[s];
		assert(pImg && pSrc);

		int numImages = pImg->getNumberOfImages();
		if (numImages != pSrc->getNumberOfImages())
		{
			setErrorString("All image and source data sets should contain the same amount of images");
			return false;
		}

		assert(s < m_numTotalPoints.size());
		m_thetas[s].resize(m_numTotalPoints[s]);
		m_betas[s].resize(m_numTotalPoints[s]);

		int pos = 0;

		for (int j = 0 ; j < numImages ; j++)
		{
			int numPoints = pImg->getNumberOfImagePoints(j);
			if (numPoints != pSrc->getNumberOfImagePoints(j))
			{
				setErrorString("All image and source data sets should have the same number of image points for all images");
				return false;
			}

			for (int k = 0 ; k < numPoints ; k++, pos++)
			{
				Vector2Dd theta = pImg->getImagePointPosition(j,k);
				Vector2Dd beta = pSrc->getImagePointPosition(j, k);
				
				theta /= m_angularScale;
				beta /= m_angularScale;
				
				m_thetas[s][pos] = Vector2Df(theta.getX(), theta.getY());
				m_betas[s][pos] = Vector2Df(beta.getX(), beta.getY());
			}			
		}
	}

	return true;
}

} // end namespace

