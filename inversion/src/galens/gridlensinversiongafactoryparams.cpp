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
#include <serut/memoryserializer.h>
#include <serut/dummyserializer.h>
#include <vector>
#include <memory>

#include "debugnew.h"

using namespace std;

namespace grale
{

GridLensInversionGAFactoryParams::GridLensInversionGAFactoryParams()
{
	m_pParams = new GridLensInversionParameters();
}

GridLensInversionGAFactoryParams::GridLensInversionGAFactoryParams(int maxgen, 
		const vector<shared_ptr<ImagesDataExtended>> &images, 
		const vector<GridSquare> &gridsquares,
		double D_d, double z_d, double massscale, 
		bool useweights,
		GridLensInversionParameters::BasisFunctionType b, bool allowNegativeValues,
		const GravitationalLens *pBaseLens,
		GridLensInversionParameters::MassSheetSearchType sheetSearchType,
		const ConfigurationParameters *pFitnessObjectParams,
		const ScaleSearchParameters &massScaleSearchParams)
{
	shared_ptr<GravitationalLens> sheetLens = GridLensInversionParameters::createDefaultSheetLens(sheetSearchType, D_d);
	m_pParams = new GridLensInversionParameters(maxgen, images, gridsquares, D_d, z_d, massscale,
			                                    useweights, b, allowNegativeValues, pBaseLens, sheetLens.get(),
												pFitnessObjectParams, massScaleSearchParams);
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

} // end namespace

