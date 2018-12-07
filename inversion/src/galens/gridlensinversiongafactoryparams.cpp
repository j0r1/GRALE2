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
#include "gridlensinversiongafactoryparams.h"
#include "imagesdataextended.h"
#include "gravitationallens.h"
#include "configurationparameters.h"
#include "plummerlens.h"
#include "squarelens.h"
#include "gausslens.h"
#include <serut/memoryserializer.h>
#include <serut/dummyserializer.h>
#include <vector>

#include "debugnew.h"

using namespace std;

namespace grale
{

GridLensInversionGAFactoryParams::GridLensInversionGAFactoryParams() : LensInversionGAFactoryParams(LensInversionGAFactoryParams::GridInversion)
{
	m_pParams = new GridLensInversionParameters();
}

GridLensInversionGAFactoryParams::GridLensInversionGAFactoryParams(int maxgen, 
		const vector<ImagesDataExtended *> &images, 
		const vector<GridSquare> &gridsquares,
		double D_d, double z_d, double massscale, 
		bool copyimages, bool useweights,
		GridLensInversionParameters::BasisFunctionType b, bool allowNegativeValues,
		const GravitationalLens *pBaseLens,
		GridLensInversionParameters::MassSheetSearchType sheetSearchType,
		const ConfigurationParameters *pFitnessObjectParams,
		bool wideSearch) 
	: LensInversionGAFactoryParams(LensInversionGAFactoryParams::GridInversion)
{
	m_pParams = new GridLensInversionParameters(maxgen, images, gridsquares, D_d, z_d, massscale, copyimages,
			                                    useweights, b, allowNegativeValues, pBaseLens, sheetSearchType,
												pFitnessObjectParams, wideSearch);
}

GridLensInversionGAFactoryParams::~GridLensInversionGAFactoryParams()
{
	delete m_pParams;
}
	
bool GridLensInversionGAFactoryParams::write(serut::SerializationInterface &si) const
{ 
	if (m_pParams->write(si))
		return true;
	setErrorString(m_pParams->getErrorString());
	return false;
}

bool GridLensInversionGAFactoryParams::read(serut::SerializationInterface &si)
{ 
	if (m_pParams->read(si))
		return true;
	setErrorString(m_pParams->getErrorString());
	return false;
}

GridLensInversionGAFactoryParams *GridLensInversionGAFactoryParams::createCopy() const
{
	auto pParams = m_pParams->createCopy();

	GridLensInversionGAFactoryParams *pRetVal = new GridLensInversionGAFactoryParams();
	delete pRetVal->m_pParams;
	pRetVal->m_pParams = pParams;
	return pRetVal;
}

vector<GridLensInversionGAFactoryParams::BasisLensInfo> GridLensInversionGAFactoryParams::getBasisLenses() const
{
	const vector<GridSquare> &squares = getGridSquares();
	const int numMasses = squares.size();
	vector<double> massWeights(numMasses);

	if (numMasses == 0)
	{
		setErrorString("No basis lenses are present");
		return vector<BasisLensInfo>(); // return empty vector
	}

	if (useMassWeights())
	{
		//sendMessage("Using mass weights");
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

	vector<BasisLensInfo> basisLenses;
	double massScale = getMassScale();
	double D_d = getD_d();
	int squareNumber = 0;

	for (auto squareIt = squares.begin() ; squareIt != squares.end() ; squareIt++, squareNumber++)
	{
		double gridSize = squareIt->getSize();
		Vector2D<double> gridCenter = squareIt->getCenter();
		bool returnValue;
		shared_ptr<GravitationalLens> pLens;
	
		if (getBasisFunctionType() == GridLensInversionParameters::PlummerBasis)
		{
			PlummerLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = shared_ptr<GravitationalLens>(new PlummerLens());
			returnValue = pLens->init(D_d, &lensParams);
		}
		else if (getBasisFunctionType() == GridLensInversionParameters::GaussBasis)
		{
			GaussLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = shared_ptr<GravitationalLens>(new GaussLens());
			returnValue = pLens->init(D_d, &lensParams);
		}
		else // Squares
		{
			SquareLensParams lensParams(massScale * massWeights[squareNumber], gridSize);

			pLens = shared_ptr<GravitationalLens>(new SquareLens());
			returnValue = pLens->init(D_d, &lensParams);
		}

		if (!returnValue)
		{
			setErrorString(std::string("Couldn't initialize basis lens: ") + pLens->getErrorString());
			return vector<BasisLensInfo>(); // return empty vector
		}

		basisLenses.push_back({ pLens, gridCenter, massScale*massWeights[squareNumber] });
	}

	return basisLenses;
}

} // end namespace

