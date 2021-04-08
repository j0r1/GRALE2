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
#include "lensinversiongafactorycommon.h"
#include "randomnumbergenerator.h"
#include "backprojectmatrix.h"
#include "vector2d.h"
#include "lensinversiongafactoryparamssingleplanecpu.h"
#include <vector>
#include <memory>

namespace grale
{

class LensInversionGAFactoryParamsSinglePlaneCPU;
class LensInversionGenome;
class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT LensInversionGAFactorySinglePlaneCPU : public LensInversionGAFactoryCommon
{
public:
	LensInversionGAFactorySinglePlaneCPU();
	~LensInversionGAFactorySinglePlaneCPU();

	bool init(const mogal::GAFactoryParams *p);

	GravitationalLens *createLens(const std::vector<float> &basisFunctionWeights,
	                                      const std::vector<float> &sheetValues,
										  float scaleFactor,
										  std::string &errStr) const override;

	bool initializeNewCalculation(const std::vector<float> &masses, const std::vector<float> &sheetValues) override;
	bool calculateMassScaleFitness(float scaleFactor, float &fitness) override;
	bool calculateTotalFitness(float scaleFactor, float *pFitnessValues) override;
private:
	void clear();
	bool localSubInit(double z_d, const std::vector<std::shared_ptr<ImagesDataExtended>> &images, 
	                  const std::vector<std::pair<std::shared_ptr<GravitationalLens>, Vector2D<double> > > &basisLenses,
					  const GravitationalLens *pBaseLens, const GravitationalLens *pSheetLens, 
					  const ConfigurationParameters *pFitnessObjectParams);

	std::unique_ptr<LensInversionGAFactoryParamsSinglePlaneCPU> m_pCurrentParams;

	std::vector<std::pair<std::shared_ptr<GravitationalLens>, Vector2D<double> > > m_basisLenses;
	std::shared_ptr<GravitationalLens> m_sheetLens;

	std::unique_ptr<DeflectionMatrix> m_pDeflectionMatrix;
	std::shared_ptr<BackProjectMatrix> m_pShortBPMatrix, m_pTotalBPMatrix;

	float m_sheetScale; // Stored when initializeNewCalculation is called
};

} // end namespace

