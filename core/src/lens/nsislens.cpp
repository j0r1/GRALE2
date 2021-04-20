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
#include "nsislens.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> NSISLensParams::createCopy() const
{
	return std::make_unique<NSISLensParams>(dispersion,corerad);
}

bool NSISLensParams::write(serut::SerializationInterface &si) const
{ 
	double arr[2]; 
	
	arr[0] = dispersion; 
	arr[1] = corerad; 
	
	if (!si.writeDoubles(arr,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool NSISLensParams::read(serut::SerializationInterface &si)
{ 
	double arr[2]; 
	
	if (!si.readDoubles(arr,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	dispersion = arr[0]; 
	corerad = arr[1]; 
	
	return true; 
}

NSISLens::NSISLens() : SymmetricLens(GravitationalLens::NSIS)
{
}

NSISLens::~NSISLens()
{
}

bool NSISLens::processParameters(const GravitationalLensParams *params)
{
	const NSISLensParams *p = dynamic_cast<const NSISLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'NSISLensParams'");
		return false;
	}
	
	double sigma = p->getVelocityDispersion();

	m_massFactor = CONST_PI*sigma*sigma*getLensDistance()/CONST_G;
	m_densFactor = sigma*sigma/(2.0*CONST_G*getLensDistance());
	m_core = p->getAngularCoreRadius();
	m_coreSquared = m_core*m_core;
	
	return true;
}

double NSISLens::getMassInside(double thetaLength) const
{
	return m_massFactor*(SQRT(thetaLength*thetaLength+m_coreSquared)-m_core);
}

double NSISLens::getProfileSurfaceMassDensity(double thetaLength) const
{	
	return m_densFactor/SQRT(thetaLength*thetaLength+m_coreSquared);
}

} // end namespace

