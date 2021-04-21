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
#include "sielens.h"
#include "constants.h"

namespace grale
{
	
std::unique_ptr<GravitationalLensParams> SIELensParams::createCopy() const
{
	return std::make_unique<SIELensParams>(V,F);
}

bool SIELensParams::write(serut::SerializationInterface &si) const
{ 
	double array[2]; 
	
	array[0] = V; 
	array[1] = F; 
	
	if (!si.writeDoubles(array,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool SIELensParams::read(serut::SerializationInterface &si)
{ 
	double array[2]; 
	
	if (!si.readDoubles(array,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	V = array[0]; 
	F = array[1]; 
	
	return true; 
}

SIELens::SIELens() : GravitationalLens(GravitationalLens::SIE)
{
}

SIELens::~SIELens()
{
}

bool SIELens::processParameters(const GravitationalLensParams *params)
{
	const SIELensParams *p = dynamic_cast<const SIELensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'SIELensParams'");
		return false;
	}
	
	V = p->getVelocityDispersion();
	F = p->getEllipticity();
	
	return true;
}

bool SIELens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *alpha) const
{
	double ff = SQRT(1.0-F*F);
	double factor = (4.0*CONST_PI*V*V)/(SPEED_C*SPEED_C)*SQRT(F)/ff;
	double x,y;

	x = factor*ASINH((ff/F)*theta.getX()/theta.getLength());
	y = factor*ASIN(ff*theta.getY()/theta.getLength());
	*alpha = Vector2D<double>(x,y);
	return true;
}

bool SIELens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double factor = (4.0*CONST_PI*V*V)/(SPEED_C*SPEED_C)*SQRT(F)/(theta.getLengthSquared()*SQRT(F*F*theta.getY()*theta.getY()+theta.getX()*theta.getX()));

	axx = factor*theta.getY()*theta.getY();
	ayy = factor*theta.getX()*theta.getX();
	axy = -factor*theta.getX()*theta.getY();

	return true;
}

double SIELens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double dens;
	
	dens = (V*V/(2.0L*CONST_G*getLensDistance()))*SQRT(F/(theta.getX()*theta.getX()+F*F*theta.getY()*theta.getY()));
	
	return dens;
}

bool SIELens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	Vector2D<double> alpha;

	getAlphaVector(theta, &alpha);

	double potential = theta.getX()*alpha.getX()+theta.getY()*alpha.getY();
	
	*pPotentialValue = (D_ds/D_s)*potential;
	return true;
}

} // end namespace

