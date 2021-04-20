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
#include "pointmasslens.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> PointmassLensParams::createCopy() const
{
	return std::make_unique<PointmassLensParams>(lensmass);
}

bool PointmassLensParams::write(serut::SerializationInterface &si) const
{ 
	if (!si.writeDouble(lensmass)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool PointmassLensParams::read(serut::SerializationInterface &si)	
{ 
	if (!si.readDouble(&lensmass)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

PointmassLens::PointmassLens() : SymmetricLens(GravitationalLens::Pointmass)
{
}

PointmassLens::~PointmassLens()
{
}

bool PointmassLens::processParameters(const GravitationalLensParams *params)
{
	const PointmassLensParams *p = dynamic_cast<const PointmassLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'PointmassLensParams'");
		return false;
	}
	
	PointmassLens::mass = p->getLensMass();
	return true;
}

double PointmassLens::getMassInside(double thetaLength) const
{
	return mass;
}

double PointmassLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	return 0;
}

} // end namespace

