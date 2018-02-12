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
#include "nfwlens.h"
#include "constants.h"

#include <iostream>

namespace grale
{

NFWLensParams::NFWLensParams()
{
}

NFWLensParams::NFWLensParams(double rho_s, double theta_s)
{
	m_densityScale3D = rho_s;
	m_angularRadiusScale = theta_s;
}

NFWLensParams::~NFWLensParams()
{
}

bool NFWLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_densityScale3D) && si.writeDouble(m_angularRadiusScale)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool NFWLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_densityScale3D) && si.readDouble(&m_angularRadiusScale)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

GravitationalLensParams *NFWLensParams::createCopy() const
{
	return new NFWLensParams(m_densityScale3D, m_angularRadiusScale);
}

NFWLens::NFWLens() : SymmetricLens(GravitationalLens::NFW)
{
}

NFWLens::~NFWLens()
{
}

bool NFWLens::processParameters(const GravitationalLensParams *p)
{
	const NFWLensParams *pParams = dynamic_cast<const NFWLensParams *>(p);

	if (!pParams)
	{
		setErrorString("Parameters are not of type 'NFWLensParams'");
		return false;
	}

	double rho_s = pParams->get3DDensityScale();
	double theta_s = pParams->getAngularRadiusScale();
	double r_s = theta_s * getLensDistance();
	
	m_angularRadiusScale = theta_s;
	m_massScale = 4.0*CONST_PI*r_s*r_s*r_s*rho_s;
	m_densityScale = 2.0*r_s*rho_s;
	m_deflectionScale = 16.0*CONST_PI*CONST_G*r_s*r_s*rho_s/(SPEED_C*SPEED_C);

	//std::cerr << "deflection scale = " << m_deflectionScale/ANGLE_ARCSEC << " arcsec" << std::endl;
	return true;
}

double NFWLens::getMassInside(double thetaLength) const
{
	double x = thetaLength/m_angularRadiusScale;

	return m_massScale*(LN(0.5*x)+F(x));
}

bool NFWLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double thetaLength = theta.getLength();
	double x = thetaLength/m_angularRadiusScale;
	
	double tx = theta.getX();
	double ty = theta.getY();
	double tx2 = tx*tx;
	double ty2 = ty*ty;
	double t2 = tx2+ty2;
	double t3 = t2*SQRT(t2);
	double t4 = t2*t2;
	double t_s = m_angularRadiusScale;

	double LNplusFx = (LN(0.5*x)+F(x));
	double DFx = DF(x);

	double Axx = (ty2-tx2)/t4*LNplusFx + tx2/t4 + tx2/(t3*t_s)*DFx;
	double Ayy = (tx2-ty2)/t4*LNplusFx + ty2/t4 + ty2/(t3*t_s)*DFx;
	double Axy = -2.0*tx*ty/t4*LNplusFx + tx*ty/t4 + tx*ty/(t3*t_s)*DFx;

	axx = m_deflectionScale*t_s*Axx;
	ayy = m_deflectionScale*t_s*Ayy;
	axy = m_deflectionScale*t_s*Axy;

	return true;
}

bool NFWLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	double thetaLength = theta.getLength();
	double x = thetaLength/m_angularRadiusScale;

	*pPotentialValue = (D_ds)/(D_s) * 0.5 * m_deflectionScale * m_angularRadiusScale * P(x);
	return true;
}

double NFWLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	double x = thetaLength/m_angularRadiusScale;

	return m_densityScale*(1.0-F(x))/(x*x-1.0); // TODO: use formula by Wright & Brainerd
}

} // end namespace

