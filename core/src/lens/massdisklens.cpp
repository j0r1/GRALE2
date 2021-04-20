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
#include "massdisklens.h"

#include "debugnew.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> MassDiskLensParams::createCopy() const
{ 
	return std::make_unique<MassDiskLensParams>(m_density, m_angularRadius); 
}

bool MassDiskLensParams::write(serut::SerializationInterface &si) const
{ 
	double x[2] = { m_density, m_angularRadius }; 
	
	if (!si.writeDoubles(x, 2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	return true; 
}

bool MassDiskLensParams::read(serut::SerializationInterface &si)
{ 
	double x[2]; 
	
	if (!si.readDoubles(x, 2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	m_density = x[0]; 
	m_angularRadius = x[1]; 
	
	return true; 
}

MassDiskLens::MassDiskLens() : GravitationalLens(GravitationalLens::MassDisk)
{
}

MassDiskLens::~MassDiskLens()
{
}

bool MassDiskLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const MassDiskLensParams *pParams = dynamic_cast<const MassDiskLensParams *>(pLensParams);

	if (!pParams)
	{
		setErrorString("Parameters are not of type 'MassDiskLensParams'");
		return false;
	}

	m_density = pParams->getDensity();
	m_angularRadius = pParams->getAngularRadius();
	m_angularRadiusSquared = m_angularRadius*m_angularRadius;
	m_factor = (4.0*CONST_PI*CONST_G*m_density*getLensDistance())/(SPEED_C*SPEED_C);

	double realRadius = m_angularRadius*getLensDistance();
	double totalMass = realRadius*realRadius*CONST_PI*m_density;
	double massScale = ABS(totalMass);
	m_scaledMass = totalMass/massScale;
	m_scaleFactor = SQRT((4.0*CONST_G*massScale)/(SPEED_C*SPEED_C*getLensDistance()));
	return true;
}

bool MassDiskLens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
	double thetaSquared = theta.getLengthSquared();
	if (thetaSquared < m_angularRadiusSquared) // inside the mass disk
		*pAlpha = theta*m_factor;
	else
	{
		Vector2D<double> scaledTheta = theta/m_scaleFactor;
		*pAlpha = (theta/scaledTheta.getLengthSquared())*m_scaledMass;
	}
	return true;
}

double MassDiskLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	if (theta.getLengthSquared() < m_angularRadiusSquared)
		return m_density;
	return 0;
}

bool MassDiskLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double thetaSquared = theta.getLengthSquared();

	if (thetaSquared < m_angularRadiusSquared) // inside the mass disk
	{
		axx = m_factor;
		ayy = m_factor;
		axy = 0;
	}
	else
	{
		Vector2D<double> scaledTheta = theta/m_scaleFactor;
		double lengthSquared = scaledTheta.getLengthSquared();
		double lengthSquaredSquared = lengthSquared*lengthSquared;
		double x = scaledTheta.getX();
		double y = scaledTheta.getY();
		double x2 = x*x;
		double y2 = y*y;

		axx = m_scaledMass*(y2-x2)/(lengthSquaredSquared);
		ayy = m_scaledMass*(x2-y2)/(lengthSquaredSquared);
		axy = m_scaledMass*(-2.0*x*y)/(lengthSquaredSquared);
	}
	return true;
}
	
bool MassDiskLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	double thetaSquared = theta.getLengthSquared();

	if (thetaSquared < m_angularRadiusSquared)
		*pPotentialValue = 0.5*m_factor*(D_ds/D_s)*thetaSquared;
	else
		*pPotentialValue = 0.5*m_factor*(D_ds/D_s)*m_angularRadiusSquared*(LN(thetaSquared/m_angularRadiusSquared)+1.0);
	return true;
}

} // end namespace

