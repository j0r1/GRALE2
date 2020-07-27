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
#include "lensinversionparameterssingleplanecpu.h"
#include <mogal/gafactory.h>

namespace grale
{

class ImagesDataExtended;
class GravitationalLens;
class ConfigurationParameters;

class GRALE_IMPORTEXPORT LensInversionGAFactoryParamsSinglePlaneCPU : public mogal::GAFactoryParams
{
public:
	LensInversionGAFactoryParamsSinglePlaneCPU();
	LensInversionGAFactoryParamsSinglePlaneCPU(const LensInversionParametersSinglePlaneCPU &params);
	~LensInversionGAFactoryParamsSinglePlaneCPU();

	int getMaximumNumberOfGenerations() const										{ return m_pParams->getMaximumNumberOfGenerations(); }
	double getD_d() const															{ return m_pParams->getD_d(); }
	double getZ_d() const															{ return m_pParams->getZ_d(); }
	double getMassScale() const														{ return m_pParams->getMassScale(); }
	const std::vector<std::shared_ptr<ImagesDataExtended>> &getImages() const		{ return m_pParams->getImages(); }
	bool allowNegativeValues() const												{ return m_pParams->allowNegativeValues(); }
	const GravitationalLens *getBaseLens() const									{ return m_pParams->getBaseLens(); }
	const GravitationalLens *getSheetLens() const									{ return m_pParams->getSheetLens(); }
	const ConfigurationParameters *getFitnessObjectParameters() const				{ return m_pParams->getFitnessObjectParameters(); }
	const ScaleSearchParameters &getMassScaleSearchParameters() const				{ return m_pParams->getMassScaleSearchParameters(); }
	const std::vector<LensInversionBasisLensInfo> &getBasisLenses() const 			{ return m_pParams->getBasisLenses(); }
	
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);

	LensInversionGAFactoryParamsSinglePlaneCPU *createCopy() const;
private:
	LensInversionParametersSinglePlaneCPU *m_pParams;
};
	
} // end namespace
