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

#include "lensfitnessj1004.h"
#include "fitnesscomponent.h"
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

#define J1004_MASK_CORNERS 1
#define J1004_MASK_POINTS 2
#define J1004_MASK_CAUSTICPENALTY 4

LensFitnessJ1004Test::LensFitnessJ1004Test()
{
	m_initialized = false;
	m_pCache = 0;
}

LensFitnessJ1004Test::~LensFitnessJ1004Test()
{
	delete m_pCache;
}

bool LensFitnessJ1004Test::init(double z_d, std::list<ImagesDataExtended *> &images, 
		                        std::list<ImagesDataExtended *> &shortImages, const ConfigurationParameters *pParams)
{
	if (m_initialized)
	{
		setErrorString("Already initialized");
		return false;
	}

	m_pointGroupsAll.clear();
	m_pointGroupsShort.clear();
	
	std::list<ImagesDataExtended *>::const_iterator it;
	int count = 0;
	
	if (images.size() == 0)
	{
		setErrorString("No available images data");
		return false;
	}

	if (images.size()%3 != 0)
	{
		setErrorString("Number of images data sets should be a multiple of three");
		return false;
	}

	delete m_pCache;
	m_pCache = new FitnessComponentCache(images.size());

	int numSources = images.size()/3;

	m_nullTriangles.resize(numSources);
	m_nullTriangleAreas.resize(numSources);
	m_fitnessMask.resize(numSources);
	m_lineSegments.resize(numSources);
	m_lineSegmentFlags.resize(numSources);
	m_lineSegmentIntersections.resize(numSources);
	m_critTriangles.resize(numSources);

	m_sourceIndices.clear();
	m_shortSourceIndices.clear();
	m_nullIndices.clear();
	m_critSources.clear();
	m_critIndices.clear();
	m_rectFlags.clear();
	m_pointGroupFlags.clear();

	m_halfNumberOfSources = numSources;
	
	for (it = images.begin() ; it != images.end() ; it++, count++)
	{
		int sourcePos = count/3;
		ImagesDataExtended *pImgDat = *it;
		
		m_totalDeflectionFlags.push_back(true);

		if (count%3 == 0) // actual images
		{
			m_sourceIndices.push_back(count);
			m_shortSourceIndices.push_back(sourcePos);

			m_pointGroupsAll.add(pImgDat);
			m_pointGroupsShort.add(pImgDat);

			shortImages.push_back(pImgDat);
			m_shortDeflectionFlags.push_back(true);
			m_shortDerivativeFlags.push_back(false);
			m_shortPotentialFlags.push_back(false);

			m_totalDerivativeFlags.push_back(false);

			if (pImgDat->getNumberOfTimeDelays() > 0)
				m_totalPotentialFlags.push_back(true);
			else
				m_totalPotentialFlags.push_back(false);

			m_totalInverseMagnificationFlags.push_back(false);

			int num;
			
			if ((num = pImgDat->getNumberOfImages()) <= 1)
			{
				setErrorString("An images data set doesn't contain multiple images");
				return false;
			}

			for (int i = 0 ; i < num ; i++)
			{
				if (pImgDat->getNumberOfImagePoints(i) <= 2)
				{
					setErrorString("Each image must contain at least three points");
					return false;
				}
			}

			int numExtraParams = pImgDat->getNumberOfExtraParameters();

			if (numExtraParams != 1)
			{
				setErrorString("A fitness mask is required for each image");
				return false;
			}

			double dummy;
			
			if (!pImgDat->getExtraParameter("0", dummy))
			{
				if (!pImgDat->getExtraParameter("0", m_fitnessMask[sourcePos]))
				{
					setErrorString("Unable to get fitness mask parameter: " + pImgDat->getErrorString());
					return false;
				}
			}
			else
				m_fitnessMask[sourcePos] = (int)dummy;

			bool useRects = (m_fitnessMask[sourcePos]&J1004_MASK_CORNERS)?true:false;
			bool usePointGroups = (m_fitnessMask[sourcePos]&J1004_MASK_POINTS)?true:false;

			m_rectFlags.push_back(useRects);
			m_pointGroupFlags.push_back(usePointGroups);

			if (pImgDat->getNumberOfTimeDelays() > 0)
				m_tdIndices.push_back(count);
		}
		else if (count%3 == 1) // null space data
		{
			m_nullIndices.push_back(count);
			m_nullWeights.push_back(1.0f);

			m_totalDerivativeFlags.push_back(false);
			m_totalPotentialFlags.push_back(false);

			m_totalInverseMagnificationFlags.push_back(false);

			int numExtraParams = pImgDat->getNumberOfExtraParameters();

			if (numExtraParams > 0)
			{
				setErrorString("No extra parameters are allowed for null space data");
				return false;
			}

			if (pImgDat->getNumberOfImages() != 1)
			{
				setErrorString("Null space data should consist of only one image");
				return false;
			}

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
		else // critical line grid
		{
			if (m_fitnessMask[sourcePos]&J1004_MASK_CAUSTICPENALTY)
			{
				assert(count-2 == sourcePos*3);
				m_critSources.push_back(count-2);
				m_critIndices.push_back(count);
			}

			m_totalDerivativeFlags.push_back(true);
			m_totalPotentialFlags.push_back(false);

			m_totalInverseMagnificationFlags.push_back(true);

			int numExtraParams = pImgDat->getNumberOfExtraParameters();

			if (numExtraParams > 0)
			{
				setErrorString("No extra parameters are allowed for critical line grid");
				return false;
			}

			if (!pImgDat->hasTriangulation())
			{
				setErrorString("Critical line data doesn't contain a triangulation");
				return false;
			}
			
			int numImages = pImgDat->getNumberOfImages();

			m_lineSegments[sourcePos].resize(numImages);
			m_lineSegmentFlags[sourcePos].resize(numImages);
			m_lineSegmentIntersections[sourcePos].resize(numImages);
			m_critTriangles[sourcePos].resize(numImages);

			for (int i = 0 ; i < numImages ; i++)
			{
				std::vector<TriangleIndices> t;

				if (!pImgDat->getTriangles(i, t))
				{
					setErrorString(std::string("Unable to obtain triangulation data from the critical line data: ") + (*it)->getErrorString());
					return false;
				}

				m_critTriangles[sourcePos][i].resize(t.size());

				int j = 0;
				std::vector<TriangleIndices>::const_iterator it;

				for (it = t.begin(); it != t.end() ; it++, j++)
				{
					TriangleIndices triangle = *it;

					int lineIndex1 = storeTriangleLine(sourcePos, i, triangle.getIndex(0), triangle.getIndex(1));
					int lineIndex2 = storeTriangleLine(sourcePos, i, triangle.getIndex(1), triangle.getIndex(2));
					int lineIndex3 = storeTriangleLine(sourcePos, i, triangle.getIndex(2), triangle.getIndex(0));
				
					m_critTriangles[sourcePos][i][j] = TriangleIndices(lineIndex1, lineIndex2, lineIndex3);
				}

				m_lineSegmentFlags[sourcePos][i].resize(m_lineSegments[sourcePos][i].size());
				m_lineSegmentIntersections[sourcePos][i].resize(m_lineSegments[sourcePos][i].size());
			}
		}
	}

	m_initialized = true;
	return true;
}

int LensFitnessJ1004Test::storeTriangleLine(int sourcePos, int imgIndex, int point1, int point2)
{
	int p1 = point1;
	int p2 = point2;

	if (p1 > p2)
	{
		p1 = point2;
		p2 = point1;
	}

	for (int i = 0 ; i < m_lineSegments[sourcePos][imgIndex].size() ; i++)
	{
		if (m_lineSegments[sourcePos][imgIndex][i].first == p1 && m_lineSegments[sourcePos][imgIndex][i].second == p2)
		return i;
	}

	int index = m_lineSegments[sourcePos][imgIndex].size();

	m_lineSegments[sourcePos][imgIndex].push_back(std::pair<int,int>(p1, p2));

	return index;
}


bool LensFitnessJ1004Test::calculateMassScaleFitness(const ProjectedImagesInterface &iface, float &fitness) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}
	
	fitness = calculateOverlapFitness_Extended(m_pointGroupsShort, iface, m_shortSourceIndices, 
			                                         m_rectFlags, m_pointGroupFlags);
	
	return true;
}

bool LensFitnessJ1004Test::calculateOverallFitness(const ProjectedImagesInterface &iface, float *fitnessvalues) const
{
	if (!isInitialized())
	{
		setErrorString("Not initialized");
		return false;
	}

	m_pCache->clear();

	fitnessvalues[0] = calculateOverlapFitness_Extended(m_pointGroupsAll, iface, m_sourceIndices, m_rectFlags,
			                                            m_pointGroupFlags, m_pCache) * m_sourceIndices.size();
	fitnessvalues[1] = calculateCausticPenaltyFitness(iface, m_critSources, m_critIndices,
			                                          m_lineSegments, m_critTriangles, m_lineSegmentFlags,
													  m_lineSegmentIntersections, m_pCache);
	fitnessvalues[2] = calculateNullFitness_ExtendedImages(iface, m_sourceIndices, m_nullIndices,
			                                               m_nullTriangles, m_nullTriangleAreas,
														   m_nullWeights,
														   m_pCache);
	fitnessvalues[3] = calculateTimeDelayFitness(iface, m_tdIndices);

	return true;
}

string LensFitnessJ1004Test::getUsage() const
{
	return R"TEXT(
This module is kept only for backwards compatibility. 

Please use the 'general' module instead.

)TEXT";
}

} // end namespace

