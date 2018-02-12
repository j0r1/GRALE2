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
#include "sersiclens.h"
#include "constants.h"
#include <gsl/gsl_sf_hyperg.h>

namespace grale
{

SersicLensParams::SersicLensParams()
{
	m_centralDensity = 0;
	m_angularScale = 0;
	m_sersicIndex = 0;
}

SersicLensParams::SersicLensParams(double centalDensity, double angularScale, double index)
{
	m_centralDensity = centalDensity;
	m_angularScale = angularScale;
	m_sersicIndex = index;
}

SersicLensParams::~SersicLensParams()
{
}

bool SersicLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_centralDensity) && si.writeDouble(m_angularScale) && si.writeDouble(m_sersicIndex)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool SersicLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_centralDensity) && si.readDouble(&m_angularScale) && si.readDouble(&m_sersicIndex)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

GravitationalLensParams *SersicLensParams::createCopy() const
{
	return new SersicLensParams(m_centralDensity, m_angularScale, m_sersicIndex);
}

SersicLens::SersicLens() : SymmetricLens(GravitationalLens::Sersic)
{
}

SersicLens::~SersicLens()
{
}

bool SersicLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const SersicLensParams *pParams = dynamic_cast<const SersicLensParams *>(pLensParams);

	if (pParams == 0)
	{
		setErrorString("Parameters are not of type 'SersicLensParams'");
		return false;
	}

	double D_d = getLensDistance();
		
	m_angularScale = pParams->getAngularScale();
	m_densityScale = pParams->getCentralDensity();
	m_index = pParams->getSersicIndex();
	m_massScale = 2.0*CONST_PI*D_d*D_d*m_densityScale*m_angularScale*m_angularScale;

	return true;
}

double SersicLens::getMassInside(double thetaLength) const
{
	double x = thetaLength/m_angularScale;

	return m_massScale*x*x*Goverx2(x, m_index);
}

double SersicLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	double x = thetaLength/m_angularScale;
	double x1overN = POW(x, 1.0/m_index);

	return m_densityScale*EXP(-x1overN);
}

bool SersicLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double thetaLength = theta.getLength();
	double x = thetaLength/m_angularScale;
	double x1overN = POW(x, 1.0/m_index);
	double Goverx2Value = Goverx2(x, m_index);
	double Ax = x*Goverx2Value;
	double DAx = EXP(-x1overN) - Goverx2Value;
	double factor = (4.0*CONST_G*m_massScale/(m_angularScale*SPEED_C*SPEED_C*getLensDistance()));
	double tx = theta.getX();
	double ty = theta.getY();
	double tx2 = tx*tx;
	double ty2 = ty*ty;
	double t2 = tx2+ty2;
	double t3 = t2*thetaLength;

	axx = factor*(Ax*ty2/t3 + tx2/(m_angularScale*t2)*DAx);
	ayy = factor*(Ax*tx2/t3 + ty2/(m_angularScale*t2)*DAx);
	axy = factor*(-Ax*tx*ty/t3 + tx*ty/(m_angularScale*t2)*DAx);

	return true;
}

double SersicLens::Goverx2(double x, double n)
{
	double x1overN = POW(x, 1.0/n);
	double M = gsl_sf_hyperg_1F1(1.0, 2.0+2.0*n, x1overN);

	return EXP(-x1overN)*(x1overN*M + 2.0*n + 1.0)/(2.0*(2.0*n+1.0));
}

} // end namespace
