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
#include "gausslens.h"
#include "constants.h"
#include <cmath>

#include "debugnew.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> GaussLensParams::createCopy() const
{
	return std::make_unique<GaussLensParams>(m_mass, m_angularWidth);
}

bool GaussLensParams::write(serut::SerializationInterface &si) const
{ 
	double array[2]; 
	
	array[0] = m_mass; 
	array[1] = m_angularWidth; 
	if (!si.writeDoubles(array,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	return true; 
}

bool GaussLensParams::read(serut::SerializationInterface &si)
{ 
	double array[2]; 
	
	if (!si.readDoubles(array,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	m_mass = array[0]; 
	m_angularWidth = array[1]; 
	return true; 
}

GaussLens::GaussLens() : SymmetricLens(GravitationalLens::Gaussian)
{
}

GaussLens::~GaussLens()
{
}

bool GaussLens::processParameters(const GravitationalLensParams *params)
{
	const GaussLensParams *p = dynamic_cast<const GaussLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'GaussLensParams'");
		return false;
	}
	
	mass = p->getMass();
	angularwidth = p->getAngularWidth();
	angularwidth2 = angularwidth*angularwidth;
	return true;
}

double GaussLens::getMassInside(double thetaLength) const
{
	double t = thetaLength/angularwidth;
	double t2 = t*t;
	
	double m = mass*(1.0L-EXP(-t2/2.0L));
	return m;
}

double GaussLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	double dens = (mass/(2.0L*CONST_PI*angularwidth2*getLensDistance()*getLensDistance()))*EXP(-(thetaLength*thetaLength)/(2.0L*angularwidth2));
	return dens;
}

} // end namespace

