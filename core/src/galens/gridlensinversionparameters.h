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

#ifndef GRALE_GRIDLENSINVERSIONPARAMETERS_H

#define GRALE_GRIDLENSINVERSIONPARAMETERS_H

#include "graleconfig.h"
#include "grid.h"
#include <serut/serializationinterface.h>
#include <errut/errorbase.h>
#include <vector>

namespace grale
{

class ImagesDataExtended;
class GravitationalLens;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT GridLensInversionParameters : public errut::ErrorBase
{
public:
	enum BasisFunctionType { PlummerBasis, SquareBasis, GaussBasis };
	enum MassSheetSearchType { NoSheet, Genome, Loop };
	
	GridLensInversionParameters();
	GridLensInversionParameters(int maxgenerations,
					 const std::vector<ImagesDataExtended *> &images, 
	                 const std::vector<GridSquare> &gridsquares, 
					 double D_d,
					 double z_d,
					 double massscale,
					 bool copyimages,
					 bool useweights = false,
					 BasisFunctionType basisFunction = PlummerBasis,
					 bool allowNegativeValues = false,
					 const GravitationalLens *pBaseLens = nullptr,
					 MassSheetSearchType sheetSearchType = NoSheet,
					 const ConfigurationParameters *pFitnessObjectParams = nullptr,
					 bool wideSearch = false
					 );
	~GridLensInversionParameters();

	int getMaximumNumberOfGenerations() const										{ return m_maxGenerations; }
	double getD_d() const															{ return m_Dd; }
	double getZ_d() const															{ return m_zd; }
	double getMassScale() const														{ return m_massScale; }
	const std::vector<ImagesDataExtended *> &getImages() const						{ return m_images; }
	const std::vector<GridSquare> &getGridSquares() const								{ return m_gridSquares; }
	bool useMassWeights() const														{ return m_useMassWeights; }
	BasisFunctionType getBasisFunctionType() const									{ return m_basisFunctionType; }
	bool allowNegativeValues() const												{ return m_allowNegative; }
	const GravitationalLens *getBaseLens() const									{ return m_pBaseLens; }
	MassSheetSearchType getMassSheetSearchType() const								{ return m_massSheetSearchType; }
	const ConfigurationParameters *getFitnessObjectParameters() const				{ return m_pParams; }
	bool useWideSearch() const														{ return m_wideSearch; }
	
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);

	GridLensInversionParameters *createCopy() const;
private:
	void zero();
	void clear();

	int m_maxGenerations;
	bool m_deleteImages;
	double m_Dd, m_massScale, m_zd;
	std::vector<ImagesDataExtended *> m_images;
	std::vector<GridSquare> m_gridSquares;
	bool m_useMassWeights;
	BasisFunctionType m_basisFunctionType;
	bool m_allowNegative;
	GravitationalLens *m_pBaseLens;
	MassSheetSearchType m_massSheetSearchType;
	ConfigurationParameters *m_pParams;
	bool m_wideSearch;
};
	
} // end namespace

#endif // GRALE_GRIDLENSINVERSIONPARAMETERS_H
