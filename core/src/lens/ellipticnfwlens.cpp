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
#include "ellipticnfwlens.h"
#include "circularlensprofile.h"
#include "constants.h"

namespace grale
{

EllipticNFWLensParams::EllipticNFWLensParams()
{
	m_densityScale3D = 0;
	m_angularRadiusScale = 0;
	m_ellipticity = 0;
}

EllipticNFWLensParams::EllipticNFWLensParams(double rho_s, double theta_s, double q)
{
	m_densityScale3D = rho_s;
	m_angularRadiusScale = theta_s;
	m_ellipticity = q;
}

EllipticNFWLensParams::~EllipticNFWLensParams()
{
}

bool EllipticNFWLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_densityScale3D) && si.writeDouble(m_angularRadiusScale) && si.writeDouble(m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool EllipticNFWLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_densityScale3D) && si.readDouble(&m_angularRadiusScale) && si.readDouble(&m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> EllipticNFWLensParams::createCopy() const
{
	return std::make_unique<EllipticNFWLensParams>(m_densityScale3D, m_angularRadiusScale, m_ellipticity);
}

class CircularNFWProfile : public CircularLensProfile
{
public:
	CircularNFWProfile(double rho_s, double theta_s, double D_d)
	{
		double r_s = D_d * theta_s;

		m_massScale = 4.0*CONST_PI*r_s*r_s*r_s*rho_s;
		m_densScale = 2.0*r_s*rho_s;
		m_derivScale = 2.0*D_d*rho_s/theta_s;

		m_angularRadiusScale = theta_s;
	}

	double getMassInside(double theta)
	{
		double x = theta/m_angularRadiusScale;

		return m_massScale * ( LN(0.5*x) + F(x) );
	}

	double getSurfaceMassDensity(double theta)
	{
		double x = theta/m_angularRadiusScale;

		return m_densScale*(1.0-F(x))/(x*x-1.0);
	}

	double getSurfaceMassDensityDerivativeOverTheta(double theta)
	{
		double x = theta/m_angularRadiusScale;
		double x2 = x*x;
		double x2min1 = x2-1.0;

		return - m_derivScale * (1.0+2.0*x2-3.0*x2*F(x))/(x2*x2min1*x2min1);
	}
private:
	double F(double x)
	{
		if (x < 1.0)
		{
			double tmp = SQRT(1.0-x*x);
			return ATANH(tmp)/tmp;
		}
		else if (x > 1.0)
		{
			double tmp = SQRT(x*x-1.0);
			return ATAN(tmp)/tmp;
		}
		return 1.0;
	}

	double m_massScale;
	double m_densScale;
	double m_derivScale;
	double m_angularRadiusScale;
};

EllipticNFWLens::EllipticNFWLens() : EllipticLens(GravitationalLens::EllipticNFW)
{
}

EllipticNFWLens::~EllipticNFWLens()
{
}

bool EllipticNFWLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const EllipticNFWLensParams *pParams = dynamic_cast<const EllipticNFWLensParams *>(pLensParams);

	if (pParams == 0)
	{
		setErrorString("Lens parameters are not of type 'EllipticNFWLensParams'");
		return false;
	}

	double q = pParams->getEllipticity();
	double rho_s = pParams->get3DDensityScale();
	double theta_s = pParams->getAngularRadiusScale();

	m_pProfile = std::make_unique<CircularNFWProfile>(rho_s, theta_s, getLensDistance());
	
	subInit(q, m_pProfile.get(), 0, 1e-5, 1024);

	return true;
}

} // end namespace
