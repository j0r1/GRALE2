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
#include "symmetriclens.h"
#include "constants.h"

namespace grale
{

SymmetricLens::SymmetricLens(GravitationalLens::LensType t) : GravitationalLens(t)
{
}

SymmetricLens::~SymmetricLens()
{
}

bool SymmetricLens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
	Vector2D<double> alphavect(0, 0);
	double mass = getMassInside(theta.getLength());

	if (mass != 0)
	{
		double massScale = ABS(mass);
		double scaledMass = mass/massScale;
		double scalefact = SQRT((4.0*CONST_G*massScale)/(SPEED_C*SPEED_C*getLensDistance()));
		Vector2D<double> scaledtheta = theta/scalefact;

		double length2 = scaledtheta.getLengthSquared();

		if (length2 != 0) // no need to divide if theta == (0, 0)
			alphavect = (theta/(scaledtheta.getLengthSquared()))*scaledMass;
	}

	*pAlpha = alphavect;
	return true;
}

bool SymmetricLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double axx0 = 0;
	double axy0 = 0;
	double ayy0 = 0;
	double thetaLength = theta.getLength();
	double massInside = getMassInside(thetaLength);

	if (massInside != 0)
	{
		double massScale = ABS(massInside);
		double scaledMass = massInside/massScale;
		double scaleFactor = SQRT(4.0*CONST_G*massScale/(SPEED_C*SPEED_C*getLensDistance()));
		Vector2D<double> scaledTheta = theta/scaleFactor;
		double t2 = scaledTheta.getLengthSquared();
		double t4 = t2*t2;
		double x = scaledTheta.getX();
		double y = scaledTheta.getY();
		double x2 = x*x;
		double y2 = y*y;
		double xy = x*y;
		double part1 = 0;
		double part3 = 0;
	
		if (t4 != 0)
		{
			part1 = scaledMass*(y2-x2)/t4;
			part3 = -scaledMass*2.0*xy/t4;
		}

		axx0 += part1;
		ayy0 += -part1;
		axy0 += part3;
	}

	double factor2 = (8.0*CONST_PI*CONST_G*getLensDistance()*getProfileSurfaceMassDensity(thetaLength))/(SPEED_C*SPEED_C);
	double thetaLengthSquared = theta.getLengthSquared();

	if (thetaLengthSquared != 0)
	{
		double x = theta.getX();
		double y = theta.getY();
		double x2 = x*x;
		double y2 = y*y;
		double xy = x*y;

		double part2 = factor2/thetaLengthSquared;

		axx0 += part2*x2;
		ayy0 += part2*y2;
		axy0 += part2*xy;
	}
	else
	{
		axx0 += factor2*0.5;
		ayy0 += factor2*0.5;
		// third part should be zero in (0,0) by symmetry (choice of axis orientation doesn't change anything about the view)
	}

	axx = axx0;
	ayy = ayy0;
	axy = axy0;

	return true;
}

SymmetricLens *SymmetricLens::cast(GravitationalLens *pLens)
{
    return dynamic_cast<SymmetricLens *>(pLens);
}

const SymmetricLens *SymmetricLens::cast(const GravitationalLens *pLens)
{
    return dynamic_cast<const SymmetricLens *>(pLens);
}

} // end namespace

