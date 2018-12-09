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

#ifndef GRALE_GRIDLENSINVERSIONGAFACTORYPARAMS_H

#define GRALE_GRIDLENSINVERSIONGAFACTORYPARAMS_H

#include "graleconfig.h"
#include "lensinversiongafactoryparams.h"
#include "gridlensinversionparameters.h"

namespace grale
{

class ImagesDataExtended;
class GravitationalLens;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT GridLensInversionGAFactoryParams : public LensInversionGAFactoryParams
{
public:
	GridLensInversionGAFactoryParams();
	GridLensInversionGAFactoryParams(int maxgenerations,
					 const std::vector<ImagesDataExtended *> &images, 
	                 const std::vector<GridSquare> &gridsquares, 
					 double D_d,
					 double z_d,
					 double massscale,
					 bool copyimages,
					 bool useweights = false,
					 GridLensInversionParameters::BasisFunctionType basisFunction = GridLensInversionParameters::PlummerBasis,
					 bool allowNegativeValues = false,
					 const GravitationalLens *pBaseLens = nullptr,
					 GridLensInversionParameters::MassSheetSearchType sheetSearchType = GridLensInversionParameters::NoSheet,
					 const ConfigurationParameters *pFitnessObjectParams = nullptr,
					 bool wideSearch = false
					 );
	~GridLensInversionGAFactoryParams();

	int getMaximumNumberOfGenerations() const										{ return m_pParams->getMaximumNumberOfGenerations(); }
	double getD_d() const															{ return m_pParams->getD_d(); }
	double getZ_d() const															{ return m_pParams->getZ_d(); }
	double getMassScale() const														{ return m_pParams->getMassScale(); }
	const std::vector<ImagesDataExtended *> &getImages() const						{ return m_pParams->getImages(); }
	bool allowNegativeValues() const												{ return m_pParams->allowNegativeValues(); }
	const GravitationalLens *getBaseLens() const									{ return m_pParams->getBaseLens(); }
	GridLensInversionParameters::MassSheetSearchType getMassSheetSearchType() const	{ return m_pParams->getMassSheetSearchType(); }
	const ConfigurationParameters *getFitnessObjectParameters() const				{ return m_pParams->getFitnessObjectParameters(); }
	bool useWideSearch() const														{ return m_pParams->useWideSearch(); }

	// TODO: for now we'll generate this from the grid, but in
	//       the future it will be stored in the constructor
	std::vector<GridLensInversionParameters::BasisLensInfo> getBasisLenses() const	{ return m_pParams->getBasisLenses(); }
	
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);

	GridLensInversionGAFactoryParams *createCopy() const;
private:
	GridLensInversionParameters *m_pParams;
};
	
} // end namespace

#endif // GRALE_GRIDLENSINVERSIONGAFACTORYPARAMS_H
