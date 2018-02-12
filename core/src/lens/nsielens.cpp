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
#include "nsielens.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{

GravitationalLensParams *NSIELensParams::createCopy() const
{
	return new NSIELensParams(m_velocityDispersion, m_ellipticity, m_angularCoreRadius);
}

bool NSIELensParams::write(serut::SerializationInterface &si) const
{ 
	double array[3]; 
	
	array[0] = m_velocityDispersion; 
	array[1] = m_ellipticity; 
	array[2] = m_angularCoreRadius; 
	
	if (!si.writeDoubles(array, 3)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool NSIELensParams::read(serut::SerializationInterface &si)
{ 
	double array[3]; 
	
	if (!si.readDoubles(array, 3)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	m_velocityDispersion = array[0]; 
	m_ellipticity = array[1]; 
	m_angularCoreRadius = array[2]; 
	
	return true; 
}

NSIELens::NSIELens() : GravitationalLens(GravitationalLens::NSIE)
{
}

NSIELens::~NSIELens()
{
}

bool NSIELens::processParameters(const GravitationalLensParams *params)
{
	const NSIELensParams *p = dynamic_cast<const NSIELensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'NSIELensParams'");
		return false;
	}
	
	V = p->getVelocityDispersion();
	F = p->getEllipticity();
	corerad = p->getAngularCoreRadius();
	
	return true;
}

bool NSIELens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *alpha0) const
{
	double factor = 4.0*CONST_PI*V*V*SQRT(F)/(SPEED_C*SPEED_C*SQRT(1.0-F*F));
	double xfrac = (theta.getX()*SQRT(1.0-F*F))/(SQRT(F*F*theta.getY()*theta.getY()+corerad*corerad+theta.getX()*theta.getX())+F*corerad);
	double yfrac = (F*theta.getY()*SQRT(1.0-F*F))/(F*SQRT(F*F*theta.getY()*theta.getY()+corerad*corerad+theta.getX()*theta.getX())+corerad);

	Vector2D<double> alpha = Vector2D<double>(factor*ATANH(xfrac), factor*ATAN(yfrac));

	*alpha0 = alpha;
	
	return true;
}


bool NSIELens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double factor = 4.0*CONST_PI*V*V*SQRT(F)/(SPEED_C*SPEED_C*SQRT(1.0-F*F));
	double x = theta.getX();
	double y = theta.getY();
	double g = SQRT(F*F*y*y+corerad*corerad+x*x);
	double denom = (F*F*y*y+corerad*corerad+2.0*F*corerad*g+F*F*corerad*corerad+x*x*F*F)*g;
	
	axx = factor*SQRT(1.0-F*F)*(F*F*y*y+corerad*corerad+F*corerad*g)/denom;
	ayy = factor*F*SQRT(1.0-F*F)*(F*corerad*corerad+F*x*x+corerad*g)/denom;
	axy = -factor*x*y*F*F*SQRT(1.0-F*F)/denom;
	
	return true;
}


double NSIELens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double dens;
		
	dens = (V*V/(2.0L*CONST_G*getLensDistance()))*SQRT(F/(theta.getX()*theta.getX()+F*F*theta.getY()*theta.getY()+corerad*corerad));
	
	return dens;
}

bool NSIELens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	Vector2D<double> alpha;

	getAlphaVector(theta, &alpha);

	double factor = 2.0*CONST_PI*V*V*corerad/(SPEED_C*SPEED_C*SQRT(F));
	double part = SQRT(F*F*theta.getY()*theta.getY()+corerad*corerad+theta.getX()*theta.getX())+corerad/F;

	double potential = theta.getX()*alpha.getX()+theta.getY()*alpha.getY() 
		         - factor*LN(part*part+(1.0-F*F)*theta.getY()*theta.getY());
	
	*pPotentialValue = (D_ds/D_s)*potential;
	return true;
}

} // end namespace

