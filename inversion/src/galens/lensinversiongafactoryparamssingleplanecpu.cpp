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
#include "lensinversiongafactoryparamssingleplanecpu.h"
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

LensInversionGAFactoryParamsSinglePlaneCPU::LensInversionGAFactoryParamsSinglePlaneCPU()
{
	m_pParams = new GridLensInversionParameters();
}

LensInversionGAFactoryParamsSinglePlaneCPU::LensInversionGAFactoryParamsSinglePlaneCPU(const GridLensInversionParameters &params)
{
	m_pParams = new GridLensInversionParameters(params);
}

LensInversionGAFactoryParamsSinglePlaneCPU::~LensInversionGAFactoryParamsSinglePlaneCPU()
{
	delete m_pParams;
}
	
bool LensInversionGAFactoryParamsSinglePlaneCPU::write(serut::SerializationInterface &si) const
{ 
	if (m_pParams->write(si))
		return true;
	setErrorString(m_pParams->getErrorString());
	return false;
}

bool LensInversionGAFactoryParamsSinglePlaneCPU::read(serut::SerializationInterface &si)
{ 
	if (m_pParams->read(si))
		return true;
	setErrorString(m_pParams->getErrorString());
	return false;
}

LensInversionGAFactoryParamsSinglePlaneCPU *LensInversionGAFactoryParamsSinglePlaneCPU::createCopy() const
{
	auto pParams = m_pParams->createCopy();

	LensInversionGAFactoryParamsSinglePlaneCPU *pRetVal = new LensInversionGAFactoryParamsSinglePlaneCPU();
	delete pRetVal->m_pParams;
	pRetVal->m_pParams = pParams;
	return pRetVal;
}

} // end namespace

