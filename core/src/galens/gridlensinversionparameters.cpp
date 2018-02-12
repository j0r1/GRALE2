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
#include "gridlensinversionparameters.h"
#include "imagesdataextended.h"
#include "gravitationallens.h"
#include "configurationparameters.h"
#include <serut/memoryserializer.h>
#include <serut/dummyserializer.h>
#include <vector>

#include "debugnew.h"

using namespace std;

namespace grale
{

GridLensInversionParameters::GridLensInversionParameters()
{
	zero();
}

GridLensInversionParameters::GridLensInversionParameters(int maxgen, 
		const vector<ImagesDataExtended *> &images, 
		const vector<GridSquare> &gridsquares,
		double D_d, double z_d, double massscale, 
		bool copyimages, bool useweights,
		BasisFunctionType b, bool allowNegativeValues,
		const GravitationalLens *pBaseLens,
		MassSheetSearchType sheetSearchType,
		const ConfigurationParameters *pFitnessObjectParams,
		bool wideSearch) 
{
	zero();

	if (copyimages)
		m_deleteImages = true;
	else
		m_deleteImages = false;
	
	m_maxGenerations = maxgen;
	m_useMassWeights = useweights; 

	m_Dd = D_d;
	m_zd = z_d;
	m_massScale = massscale; 
	
	for (auto it = images.begin() ; it != images.end() ; it++)
	{
		ImagesDataExtended *img = (*it);
		ImagesDataExtended *img2;

		if (copyimages)
			img2 = new ImagesDataExtended(*img);
		else
			img2 = img;
			
		m_images.push_back(img2);
	}

	m_gridSquares = gridsquares;
	m_basisFunctionType = b;
	m_allowNegative = allowNegativeValues;

	m_pBaseLens = nullptr;
	if (pBaseLens)
		m_pBaseLens = pBaseLens->createCopy();

	m_massSheetSearchType = sheetSearchType;

	if (pFitnessObjectParams)
		m_pParams = new ConfigurationParameters(*pFitnessObjectParams);

	m_wideSearch = wideSearch;
}

GridLensInversionParameters::~GridLensInversionParameters()
{
	clear();
}

#define BASISFUNCTION_PLUMMER 1
#define BASISFUNCTION_SQUARE  2
#define BASISFUNCTION_GAUSS   3

#define SHEETSEARCHTYPE_NONE  	0
#define SHEETSEARCHTYPE_GENOME	1
#define SHEETSEARCHTYPE_LOOP	2

bool GridLensInversionParameters::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32(m_maxGenerations))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeDouble(m_Dd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeDouble(m_zd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.writeDouble(m_massScale))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t sizes[2];

	sizes[0] = m_images.size();
	sizes[1] = m_gridSquares.size();

	if (!si.writeInt32s(sizes,2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (auto pImg : m_images)
	{
		assert(pImg);

		if (!pImg->write(si))
		{
			setErrorString(pImg->getErrorString());
			return false;
		}
	}

	for (auto s : m_gridSquares)
	{
		double params[3];

		params[0] = s.getSize();
		params[1] = s.getCenter().getX();
		params[2] = s.getCenter().getY();
	
		if (!si.writeDoubles(params,3))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}

	int32_t w = (m_useMassWeights)?1:0;
	if (!si.writeInt32(w))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t b = 0;

	switch(m_basisFunctionType)
	{
	case PlummerBasis:
		b = BASISFUNCTION_PLUMMER;
		break;
	case SquareBasis:
		b = BASISFUNCTION_SQUARE;
		break;
	case GaussBasis:
		b = BASISFUNCTION_GAUSS;
		break;
	default:
		setErrorString("Detected invalid basis function setting");
		return false;
	}

	if (!si.writeInt32(b))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t n = (m_allowNegative)?1:0;

	if (!si.writeInt32(n))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t val = (m_pBaseLens != nullptr)?1:0;
	
	if (!si.writeInt32(val))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (m_pBaseLens)
	{
		if (!m_pBaseLens->write(si))
		{
			setErrorString(m_pBaseLens->getErrorString());
			return false;
		}
	}

	switch(m_massSheetSearchType)
	{
	case Loop:
		val = SHEETSEARCHTYPE_LOOP;
		break;
	case Genome:
		val = SHEETSEARCHTYPE_GENOME;
		break;
	case NoSheet:
	default:
		val = SHEETSEARCHTYPE_NONE;
	}

	if (!si.writeInt32(val))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	val = (m_pParams == nullptr)?0:1;
	if (!si.writeInt32(val)) // write marker for additional parameters
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (m_pParams)
	{
		if (!m_pParams->write(si))
		{
			setErrorString(m_pParams->getErrorString());
			return false;
		}
	}

	int32_t searchType = (m_wideSearch)?1:0;
	if (!si.writeInt32(searchType))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	return true;
}

void GridLensInversionParameters::clear()
{
	if (m_deleteImages)
	{
		for (auto pImg : m_images)
			delete pImg;
	}
	m_images.clear();
	m_gridSquares.clear();

	delete m_pBaseLens;
	delete m_pParams;

	zero();
}

void GridLensInversionParameters::zero()
{
	m_maxGenerations = 0;
	m_deleteImages = false;
	m_Dd = 0;
	m_massScale = 0;
	m_zd = 0;
	m_useMassWeights = false;
	m_basisFunctionType = PlummerBasis;
	m_allowNegative = false;
	m_pBaseLens = nullptr;
	m_massSheetSearchType = NoSheet;
	m_pParams = nullptr;
}

bool GridLensInversionParameters::read(serut::SerializationInterface &si)
{
	clear();

	m_deleteImages = true;
	
	int32_t maxgen;

	if (!si.readInt32(&maxgen))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	m_maxGenerations = (int)maxgen;
	
	if (!si.readDouble(&m_Dd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readDouble(&m_zd))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (!si.readDouble(&m_massScale))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t sizes[2];

	if (!si.readInt32s(sizes,2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	int32_t numimages = sizes[0];
	int32_t numsquares = sizes[1];
	
	for (int32_t i = 0 ; i < numimages ; i++)
	{
		ImagesDataExtended *img = new ImagesDataExtended();
		if (!img->read(si))
		{
			setErrorString(img->getErrorString());
			delete img;
			return false;
		}
		m_images.push_back(img);
	}

	for (int32_t i = 0 ; i < numsquares ; i++)
	{
		double info[3];

		if (!si.readDoubles(info,3))
		{
			setErrorString(si.getErrorString());
			return false;
		}
		m_gridSquares.push_back(GridSquare(Vector2D<double>(info[1],info[2]),info[0]));
	}
	
	int32_t w;

	if (!si.readInt32(&w))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (w == 0)
		m_useMassWeights = false;
	else
		m_useMassWeights = true;

	int32_t b;

	if (!si.readInt32(&b))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	switch(b)
	{
	case BASISFUNCTION_PLUMMER:
		m_basisFunctionType = PlummerBasis;
		break;
	case BASISFUNCTION_SQUARE:
		m_basisFunctionType = SquareBasis;
		break;
	case BASISFUNCTION_GAUSS:
		m_basisFunctionType = GaussBasis;
		break;
	default:
		setErrorString("Read invalid basis function identifier");
		return false;
	}

	int32_t n;

	if (!si.readInt32(&n))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_allowNegative = (n == 0)?false:true;

	int32_t val;

	if (!si.readInt32(&val))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (val != 0)
	{
		string errStr;

		if (!GravitationalLens::read(si, &m_pBaseLens, errStr))
		{
			setErrorString(errStr);
			return false;
		}
	}

	if (!si.readInt32(&val))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	switch(val)
	{
	case SHEETSEARCHTYPE_NONE:
		m_massSheetSearchType = NoSheet;
		break;
	case SHEETSEARCHTYPE_GENOME:
		m_massSheetSearchType = Genome;
		break;
	case SHEETSEARCHTYPE_LOOP:
		m_massSheetSearchType = Loop;
		break;
	default:
		setErrorString("Invalid mass-sheet search type");
		return false;
	}

	if (!si.readInt32(&val))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	if (val != 0)
	{
		m_pParams = new ConfigurationParameters();
		if (!m_pParams->read(si))
		{
			setErrorString(m_pParams->getErrorString());
			delete m_pParams;
			m_pParams = nullptr;
			return false;
		}
	}

	int32_t searchType;
	if (!si.readInt32(&searchType))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	if (searchType == 0)
		m_wideSearch = false;
	else
		m_wideSearch = true;

	return true;
}
	
GridLensInversionParameters *GridLensInversionParameters::createCopy() const
{
	return new GridLensInversionParameters(m_maxGenerations, m_images, m_gridSquares, m_Dd, m_zd, m_massScale,
			                                    true, m_useMassWeights, m_basisFunctionType, m_allowNegative,
											    m_pBaseLens, m_massSheetSearchType, m_pParams, m_wideSearch);
}

} // end namespace

