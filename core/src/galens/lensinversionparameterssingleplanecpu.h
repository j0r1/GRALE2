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
#include "grid.h"
#include <serut/serializationinterface.h>
#include <errut/errorbase.h>
#include <vector>
#include <memory>

namespace grale
{

class ImagesDataExtended;
class GravitationalLens;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT ScaleSearchParameters : public errut::ErrorBase
{
public:
	ScaleSearchParameters(float startFactor, float stopFactor, int numIt, int firstItSteps, int subseqItSteps);
	ScaleSearchParameters(bool wideSearch);
	ScaleSearchParameters(); // No search
	~ScaleSearchParameters();

	bool operator==(const ScaleSearchParameters &src) const;

	float getStartFactor() const 									{ return m_startFactor; }
	float getStopFactor() const 									{ return m_stopFactor; }
	int getNumberOfIterations() const 								{ return m_numIterations; }
	int getStepsOnFirstIteration() const 							{ return m_firstIterationsSteps; }
	int getStepsOnSubsequentIterations() const 						{ return m_subseqIterationSteps; }

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::string toString() const;
private:
	float m_startFactor, m_stopFactor;
	int m_firstIterationsSteps, m_subseqIterationSteps;
	int m_numIterations;
};

class GRALE_IMPORTEXPORT LensInversionParametersSinglePlaneCPU : public errut::ErrorBase
{
public:
	enum BasisFunctionType { PlummerBasis, SquareBasis, GaussBasis };
	enum MassSheetSearchType { NoSheet, Genome };
	
	class BasisLensInfo
	{
	public:
		BasisLensInfo(std::shared_ptr<GravitationalLens> &lens, Vector2Dd center, double relevantLensingMass)
			: m_pLens(lens), m_center(center), m_relevantLensingMass(relevantLensingMass)
		{
		}

		BasisLensInfo(const BasisLensInfo &src)
		{
			copyFrom(src);
		}

		BasisLensInfo &operator=(const BasisLensInfo &src)
		{
			copyFrom(src);
			return *this;
		}

		std::shared_ptr<GravitationalLens> m_pLens;
		Vector2Dd m_center;
		double m_relevantLensingMass;

	private:
		void copyFrom(const BasisLensInfo &src)
		{
			m_pLens = src.m_pLens;
			m_center = src.m_center;
			m_relevantLensingMass = src.m_relevantLensingMass;
		}
	};

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
	LensInversionParametersSinglePlaneCPU(int maxGenerations,
			const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
			const std::vector<BasisLensInfo> &basisLenses,
			double D_d,
			double z_d,
			double massScale,
			bool allowNegativeValues = false,
			const GravitationalLens *pBaseLens = nullptr,
			const GravitationalLens *pSheetLens = nullptr,
			const ConfigurationParameters *pFitnessObjectParams = nullptr,
			const ScaleSearchParameters &massScaleSearchParams = ScaleSearchParameters(false));

	LensInversionParametersSinglePlaneCPU(int maxgenerations,
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
					 const ScaleSearchParameters &massScaleSearchParams = ScaleSearchParameters(false)
					 );

	LensInversionParametersSinglePlaneCPU(const LensInversionParametersSinglePlaneCPU &src)				{ copyFrom(src); }
	~LensInversionParametersSinglePlaneCPU();

	LensInversionParametersSinglePlaneCPU &operator=(const LensInversionParametersSinglePlaneCPU &src)	{ copyFrom(src); return *this;}
	int getMaximumNumberOfGenerations() const										{ return m_maxGenerations; }
	double getD_d() const															{ return m_Dd; }
	double getZ_d() const															{ return m_zd; }
	double getMassScale() const														{ return m_massScale; }
	const std::vector<std::shared_ptr<ImagesDataExtended>> &getImages() const		{ return m_images; }
	bool allowNegativeValues() const												{ return m_allowNegative; }

	const GravitationalLens *getBaseLens() const									{ return m_pBaseLens.get(); }
	const GravitationalLens *getSheetLens() const									{ return m_pSheetLens.get(); }
	const ConfigurationParameters *getFitnessObjectParameters() const				{ return m_pParams.get(); }
	const ScaleSearchParameters &getMassScaleSearchParameters() const				{ return m_scaleSearchParams; }
	
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);

	const std::vector<BasisLensInfo> &getBasisLenses() const						{ return m_basisLenses; }

	LensInversionParametersSinglePlaneCPU *createCopy() const;

	static std::shared_ptr<GravitationalLens> createDefaultSheetLens(MassSheetSearchType t, double Dd);
private:
	void copyFrom(const LensInversionParametersSinglePlaneCPU &src);
	void commonConstructor(int maxGenerations,
			const std::vector<std::shared_ptr<ImagesDataExtended>> &images,
			double D_d,
			double z_d,
			double massScale,
			bool allowNegativeValues,
			const GravitationalLens *pBaseLens,
			const GravitationalLens *pSheetLens,
			const ConfigurationParameters *pFitnessObjectParams,
			const ScaleSearchParameters &massScaleSearchParams);

	void buildBasisLenses(const std::vector<GridSquare> &squares, BasisFunctionType basisFunctionType, bool useMassWeights);
	void zero();
	void clear();

	int m_maxGenerations;
	double m_Dd, m_massScale, m_zd;
	std::vector<std::shared_ptr<ImagesDataExtended>> m_images;
	bool m_allowNegative;
	std::shared_ptr<GravitationalLens> m_pBaseLens;
	std::shared_ptr<GravitationalLens> m_pSheetLens;
	std::shared_ptr<ConfigurationParameters> m_pParams;
	ScaleSearchParameters m_scaleSearchParams;

	std::vector<BasisLensInfo> m_basisLenses;
};
	
} // end namespace
