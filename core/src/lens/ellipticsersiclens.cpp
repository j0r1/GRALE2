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
#include "ellipticsersiclens.h"
#include "constants.h"
#include "circularlensprofile.h"
#include <gsl/gsl_sf_hyperg.h>

namespace grale
{

EllipticSersicLensParams::EllipticSersicLensParams()
{
	m_centralDensity = 0;
	m_angularScale = 0;
	m_sersicIndex = 0;
	m_ellipticity = 0;
}

EllipticSersicLensParams::EllipticSersicLensParams(double centalDensity, double angularScale, double index, double q)
{
	m_centralDensity = centalDensity;
	m_angularScale = angularScale;
	m_sersicIndex = index;
	m_ellipticity = q;
}

EllipticSersicLensParams::~EllipticSersicLensParams()
{
}

bool EllipticSersicLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_centralDensity) && si.writeDouble(m_angularScale) && si.writeDouble(m_sersicIndex) && 
	      si.writeDouble(m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool EllipticSersicLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_centralDensity) && si.readDouble(&m_angularScale) && si.readDouble(&m_sersicIndex) &&
	      si.readDouble(&m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> EllipticSersicLensParams::createCopy() const
{
	return std::make_unique<EllipticSersicLensParams>(m_centralDensity, m_angularScale, m_sersicIndex, m_ellipticity);
}

class CircularSersicProfile : public CircularLensProfile
{
public:
	CircularSersicProfile(double densScale, double angScale, double n, double D_d)
	{
		m_index = n;
		m_densityScale = densScale;
		m_angularScale = angScale;
		m_massScale = 2.0*CONST_PI*D_d*D_d*m_densityScale*m_angularScale*m_angularScale;
		m_derivScale = m_densityScale/(m_angularScale*m_angularScale);
	}

	double getMassInside(double theta) const
	{
		double x = theta/m_angularScale;

		return m_massScale*G(x, m_index);
	}

	double getSurfaceMassDensity(double theta) const
	{
		double x = theta/m_angularScale;
		double x1overN = POW(x, 1.0/m_index);

		return m_densityScale*EXP(-x1overN);
	}

	double getSurfaceMassDensityDerivativeOverTheta(double theta) const
	{
		double x = theta/m_angularScale;
		double x1overN = POW(x, 1.0/m_index);

		return -m_derivScale*x1overN/(x*x)*EXP(-x1overN)/m_index;
	}

	static double G(double x, double n)
	{
		double x1overN = POW(x, 1.0/n);
		double M = gsl_sf_hyperg_1F1(1.0, 2.0+2.0*n, x1overN);

		return x*x*EXP(-x1overN)*(x1overN*M + 2.0*n + 1.0)/(2.0*(2.0*n+1.0));
	}
private:
	double m_index;
	double m_densityScale;
	double m_angularScale;
	double m_massScale;
	double m_derivScale;
};

EllipticSersicLens::EllipticSersicLens() : EllipticLens(EllipticSersic)
{
}

EllipticSersicLens::~EllipticSersicLens()
{
}

bool EllipticSersicLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const EllipticSersicLensParams *pParams = dynamic_cast<const EllipticSersicLensParams *>(pLensParams);

	if (pParams == 0)
	{
		setErrorString("Lens parameters are not of type 'EllipticSersicLensParams'");
		return false;
	}

	double D_d = getLensDistance();
	double q = pParams->getEllipticity();	
	double angularScale = pParams->getAngularScale();
	double densityScale = pParams->getCentralDensity();
	double n = pParams->getSersicIndex();

	m_pProfile = std::make_unique<CircularSersicProfile>(densityScale, angularScale, n, D_d);

	subInit(q, m_pProfile.get(), 0, 1e-5, 1024);

	return true;
}

} // end namespace
