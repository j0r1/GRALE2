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

#pragma once

#include "graleconfig.h"
#include "lensinversionparametersbase.h"
#include "grid.h"
#include "lensinversionbasislensinfo.h"
#include "scalesearchparameters.h"
#include <vector>
#include <memory>

namespace grale
{

class ImagesDataExtended;
class GravitationalLens;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT LensInversionParametersSinglePlaneCPU : public LensInversionParametersBase
{
public:
	enum BasisFunctionType { PlummerBasis, SquareBasis, GaussBasis };
	enum MassSheetSearchType { NoSheet, Genome };
	
	LensInversionParametersSinglePlaneCPU();

	// The idea is to supply a vector of basis functions, which basically 
	// corresponds to the contribution when the value in the genome is 1.
	// Each basis function also has a relevant lensing mass, which can be
	// the entire mass of the basis function, or the mass in the lensing
	// are if more appropriate (e.g. for a mass sheet, or for a SIS lens)
	// For a given genome (basis function weights), the weights will first
	// be rescaled so that the total mass corresponds to the mass scale, and
	// then a search is done for another scale factor that provides the best
	// fitness value


	LensInversionParametersSinglePlaneCPU(
			const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
			const std::vector<LensInversionBasisLensInfo> &basisLenses,
			double D_d,
			double z_d,
			double massScale,
			bool allowNegativeValues = false,
			const GravitationalLens *pBaseLens = nullptr,
			const GravitationalLens *pSheetLens = nullptr,
			const ConfigurationParameters *pFitnessObjectParams = nullptr,
			const ScaleSearchParameters &massScaleSearchParams = ScaleSearchParameters(false),
			bool randomizeImagePositions = false,
			uint64_t initialUncertSeed = 0
			);

	LensInversionParametersSinglePlaneCPU(
					 const std::vector<std::shared_ptr<ImagesDataExtended>> &images, 
	                 const std::vector<GridSquare> &gridsquares, 
					 double D_d,
					 double z_d,
					 double massscale,
					 bool useweights = false,
					 BasisFunctionType basisFunction = PlummerBasis,
					 bool allowNegativeValues = false,
					 const GravitationalLens *pBaseLens = nullptr,
					 const GravitationalLens *pSheetLens = nullptr,
					 const ConfigurationParameters *pFitnessObjectParams = nullptr,
					 const ScaleSearchParameters &massScaleSearchParams = ScaleSearchParameters(false),
					 bool randomizeImagePositions = false,
					 uint64_t initialUncertSeed = 0
					 );

	LensInversionParametersSinglePlaneCPU(const LensInversionParametersSinglePlaneCPU &src)				{ copyFrom(src); }
	~LensInversionParametersSinglePlaneCPU();

	LensInversionParametersSinglePlaneCPU &operator=(const LensInversionParametersSinglePlaneCPU &src)	{ copyFrom(src); return *this;}
	double getD_d() const															{ return m_Dd; }
	double getZ_d() const															{ return m_zd; }
	double getMassScale() const														{ return m_massScale; }
	const std::vector<std::shared_ptr<ImagesDataExtended>> &getImages() const		{ return m_images; }
	bool allowNegativeValues() const												{ return m_allowNegative; }

	const GravitationalLens *getBaseLens() const									{ return m_pBaseLens.get(); }
	const GravitationalLens *getSheetLens() const									{ return m_pSheetLens.get(); }
	const ConfigurationParameters *getFitnessObjectParameters() const				{ return m_pParams.get(); }
	const ScaleSearchParameters &getMassScaleSearchParameters() const				{ return m_scaleSearchParams; }
	bool getRandomizeInputPositions() const											{ return m_randomizeImagePositions; }
	uint64_t getInitialPositionUncertaintySeed() const								{ return m_initialUncertSeed; }

	bool write(serut::SerializationInterface &si) const override;
	bool read(serut::SerializationInterface &si) override;

	const std::vector<LensInversionBasisLensInfo> &getBasisLenses() const			{ return m_basisLenses; }

	std::unique_ptr<LensInversionParametersSinglePlaneCPU> createCopy() const;

	static std::shared_ptr<GravitationalLens> createDefaultSheetLens(MassSheetSearchType t, double Dd);
private:
	void copyFrom(const LensInversionParametersSinglePlaneCPU &src);
	void commonConstructor(
			const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
			double D_d,
			double z_d,
			double massScale,
			bool allowNegativeValues,
			const GravitationalLens *pBaseLens,
			const GravitationalLens *pSheetLens,
			const ConfigurationParameters *pFitnessObjectParams,
			const ScaleSearchParameters &massScaleSearchParams,
			bool randomizeImagePositions,
			uint64_t initialUncertSeed
			);

	void buildBasisLenses(const std::vector<GridSquare> &squares, BasisFunctionType basisFunctionType, bool useMassWeights);
	void zero();
	void clear();
	void printBasisLenses();

	double m_Dd, m_massScale, m_zd;
	std::vector<std::shared_ptr<ImagesDataExtended>> m_images;
	bool m_allowNegative;
	std::shared_ptr<GravitationalLens> m_pBaseLens;
	std::shared_ptr<GravitationalLens> m_pSheetLens;
	std::shared_ptr<ConfigurationParameters> m_pParams;
	ScaleSearchParameters m_scaleSearchParams;
	bool m_randomizeImagePositions;
	uint64_t m_initialUncertSeed;

	std::vector<LensInversionBasisLensInfo> m_basisLenses;

};
	
} // end namespace
