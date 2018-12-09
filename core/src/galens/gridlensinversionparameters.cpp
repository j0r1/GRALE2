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
#include "plummerlens.h"
#include "squarelens.h"
#include "gausslens.h"
#include <assert.h>
#include <vector>

#include "debugnew.h"

using namespace std;

namespace grale
{

GridLensInversionParameters::GridLensInversionParameters()
{
	zero();
}

void GridLensInversionParameters::commonConstructor(int maxGenerations,
			const std::vector<ImagesDataExtended *> &images,
			double D_d,
			double z_d,
			double massScale,
			bool copyImages,
			bool allowNegativeValues,
			const GravitationalLens *pBaseLens,
			MassSheetSearchType sheetSearchType,
			const ConfigurationParameters *pFitnessObjectParams,
			bool wideSearch)
{
	zero();

	if (copyImages)
		m_deleteImages = true;
	else
		m_deleteImages = false;
	
	m_maxGenerations = maxGenerations;
	m_allowNegative = allowNegativeValues;

	m_Dd = D_d;
	m_zd = z_d;
	m_massScale = massScale; 
	
	for (auto it = images.begin() ; it != images.end() ; it++)
	{
		ImagesDataExtended *img = (*it);
		ImagesDataExtended *img2;

		if (copyImages)
			img2 = new ImagesDataExtended(*img);
		else
			img2 = img;
			
		m_images.push_back(img2);
	}

	m_pBaseLens = nullptr;
	if (pBaseLens)
		m_pBaseLens = pBaseLens->createCopy();

	m_massSheetSearchType = sheetSearchType;

	if (pFitnessObjectParams)
		m_pParams = new ConfigurationParameters(*pFitnessObjectParams);

	m_wideSearch = wideSearch;
}

GridLensInversionParameters::GridLensInversionParameters(int maxGenerations,
			const std::vector<ImagesDataExtended *> &images,
			const std::vector<BasisLensInfo> &basisLenses,
			double D_d,
			double z_d,
			double massScale,
			bool copyImages,
			bool allowNegativeValues,
			const GravitationalLens *pBaseLens,
			MassSheetSearchType sheetSearchType,
			const ConfigurationParameters *pFitnessObjectParams,
			bool wideSearch)
{
	commonConstructor(maxGenerations, images, D_d, z_d, massScale, copyImages, allowNegativeValues, pBaseLens,
			          sheetSearchType, pFitnessObjectParams, wideSearch);

	m_basisLenses = basisLenses;
}

GridLensInversionParameters::GridLensInversionParameters(int maxGenerations, 
		const vector<ImagesDataExtended *> &images, 
		const vector<GridSquare> &gridsquares,
		double D_d, double z_d, double massScale, 
		bool copyImages, bool useweights,
		BasisFunctionType b, bool allowNegativeValues,
		const GravitationalLens *pBaseLens,
		MassSheetSearchType sheetSearchType,
		const ConfigurationParameters *pFitnessObjectParams,
		bool wideSearch) 
{
	commonConstructor(maxGenerations, images, D_d, z_d, massScale, copyImages, allowNegativeValues, pBaseLens,
			          sheetSearchType, pFitnessObjectParams, wideSearch);

	buildBasisLenses(gridsquares, b, useweights);
}

GridLensInversionParameters::~GridLensInversionParameters()
{
	clear();
}

#define SHEETSEARCHTYPE_NONE  	0
#define SHEETSEARCHTYPE_GENOME	1

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
	sizes[1] = m_basisLenses.size();

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

	for (const auto &bl : m_basisLenses)
	{
		assert(bl.m_pLens.get());
		if (!bl.m_pLens->write(si))
		{
			setErrorString("Unable to write a basis lens model: " + bl.m_pLens->getErrorString());
			return false;
		}

		const double centerAndMass[3] = { bl.m_center.getX(), bl.m_center.getY(), bl.m_relevantLensingMass };
		if (!si.writeDoubles(centerAndMass, 3))
		{
			setErrorString("Unable to write basis lens model position or mass: " + si.getErrorString());
			return false;
		}
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
	m_basisLenses.clear();

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
	int32_t numBasisFunctions = sizes[1];
	
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

	for (int32_t i = 0 ; i < numBasisFunctions ; i++)
	{
		string errStr;
		GravitationalLens *pBasisFunction = nullptr;
		if (!GravitationalLens::read(si, &pBasisFunction, errStr))
		{
			setErrorString("Unable to read basis lens model: " + errStr);
			return false;
		}
		shared_ptr<GravitationalLens> basisLens(pBasisFunction);

		double centerAndMass[3];
		if (!si.readDoubles(centerAndMass, 3))
		{
			setErrorString("Unable to read center or mass of a basis lens model:" + si.getErrorString());
			return false;
		}

		m_basisLenses.push_back( { basisLens, { centerAndMass[0], centerAndMass[1] }, centerAndMass[2] });
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
	vector<BasisLensInfo> copiedBasisLenses;

	for (const auto bl : m_basisLenses)
	{
		assert(bl.m_pLens);
		shared_ptr<GravitationalLens> newLens(bl.m_pLens->createCopy());
		if (!newLens.get())
		{
			setErrorString("Unable to create a copy from one of the basis lenses: " + bl.m_pLens->getErrorString());
			return nullptr;
		}

		copiedBasisLenses.push_back( { newLens, bl.m_center, bl.m_relevantLensingMass } );
	}

	return new GridLensInversionParameters(m_maxGenerations, m_images, copiedBasisLenses,
			                               m_Dd, m_zd, m_massScale, true, m_allowNegative,
										   m_pBaseLens, m_massSheetSearchType, m_pParams,
										   m_wideSearch);
}


void GridLensInversionParameters::buildBasisLenses(const vector<GridSquare> &squares,
                                                   BasisFunctionType basisFunctionType,
												   bool useMassWeights)
{
	const int numMasses = squares.size();
	vector<double> massWeights(numMasses);

	if (numMasses == 0)
		return;

	if (useMassWeights)
	{
		cerr << "Using mass weights" << endl;
		
		double sum = 0;

		auto it = squares.begin();
		for (it = squares.begin(); it != squares.end() ; it++)
		{
			double w = (*it).getSize();
			double w2 = (w*w);

			sum += w2;
		}
		
		// on a uniform grid, we'd like each weight to be one
		it = squares.begin();
		for (int i = 0 ; it != squares.end() ; i++, it++)
		{
			double w = (*it).getSize();
			double w2 = (w*w);

			massWeights[i] = ((((double)numMasses)*w2)/sum);
		}
	}
	else
	{
		int i = 0;

		for (auto it = squares.begin(); it != squares.end() ; it++, i++)
			massWeights[i] = 1;
	}

	// Create the basis functions

	double massScale = getMassScale();
	double D_d = getD_d();
	int squareNumber = 0;
	vector<BasisLensInfo> basisLenses;

	for (auto squareIt = squares.begin() ; squareIt != squares.end() ; squareIt++, squareNumber++)
	{
		double gridSize = squareIt->getSize();
		Vector2D<double> gridCenter = squareIt->getCenter();
		bool returnValue;
		shared_ptr<GravitationalLens> pLens;
	
		if (basisFunctionType == GridLensInversionParameters::PlummerBasis)
		{
			PlummerLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = shared_ptr<GravitationalLens>(new PlummerLens());
			returnValue = pLens->init(D_d, &lensParams);
		}
		else if (basisFunctionType == GridLensInversionParameters::GaussBasis)
		{
			GaussLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = shared_ptr<GravitationalLens>(new GaussLens());
			returnValue = pLens->init(D_d, &lensParams);
		}
		else // Squares
		{
			assert(basisFunctionType == GridLensInversionParameters::SquareBasis);

			SquareLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = shared_ptr<GravitationalLens>(new SquareLens());
			returnValue = pLens->init(D_d, &lensParams);
		}

		if (!returnValue)
		{
			cerr << "Couldn't initialize basis lens: " << pLens->getErrorString() << endl;
			return;
		}

		basisLenses.push_back({ pLens, gridCenter, massScale*massWeights[squareNumber] });
	}

	m_basisLenses = basisLenses;
}

} // end namespace

