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
#include "pimdlens.h"
#include "constants.h"

namespace grale
{

PIMDLensParams::PIMDLensParams() 
	: m_sigma0(0),
	  m_coreRadius(0),
	  m_scaleRadius(0)
{
}

PIMDLensParams::PIMDLensParams(double sigma0, double coreRadius, double scaleRadius)
	: m_sigma0(sigma0),
	  m_coreRadius(coreRadius),
	  m_scaleRadius(scaleRadius)
{
}

PIMDLensParams::~PIMDLensParams()
{
}

bool PIMDLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_sigma0) &&
	      si.writeDouble(m_coreRadius) &&
		  si.writeDouble(m_scaleRadius)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool PIMDLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_sigma0) &&
		  si.readDouble(&m_coreRadius) &&
		  si.readDouble(&m_scaleRadius)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> PIMDLensParams::createCopy() const
{
	return std::make_unique<PIMDLensParams>(m_sigma0, m_coreRadius, m_scaleRadius);
}


PIMDLens::PIMDLens() : SymmetricLens(GravitationalLens::PIMD)
{
}

PIMDLens::~PIMDLens()
{
}

bool PIMDLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const PIMDLensParams *pParams = dynamic_cast<const PIMDLensParams *>(pLensParams);

	if (pParams == 0)
	{
		setErrorString("Parameters are not of type 'PIMDLensParams'");
		return false;
	}

	m_sigma0 = pParams->getCentralDensity();

	m_coreRadius = pParams->getCoreRadius();
	m_scaleRadius = pParams->getScaleRadius();
	if (m_coreRadius < 0 || m_scaleRadius < m_coreRadius)
	{
		setErrorString("Invalid radii");
		return false;
	}

	double Dd = getLensDistance();
	m_massFactor = 2.0*CONST_PI*m_sigma0*m_scaleRadius*m_coreRadius/(m_scaleRadius-m_coreRadius) *Dd*Dd;
	m_densFactor = m_sigma0*m_scaleRadius*m_coreRadius/(m_scaleRadius-m_coreRadius);
	m_a2 = m_coreRadius*m_coreRadius;
	m_s2 = m_scaleRadius*m_scaleRadius;

	return true;
}

/* This code also works getAlphaVector:

	double Dd = getLensDistance();
	double Q = 8.0*CONST_PI*CONST_G/(SPEED_C*SPEED_C) * m_densFactor *Dd;

	auto f = [](double t2, double a_s) {
		return 1.0/(a_s + std::sqrt(a_s*a_s + t2));
	};

	double t2 = theta.getLengthSquared();
	double common = Q*(f(t2, m_coreRadius) - f(t2, m_scaleRadius));
	*pAlpha = Vector2Dd(common*theta.getX(), common*theta.getY());

   And this works for getAlphaVectorSecondDerivatives

	double Dd = getLensDistance();
	double Q = 8.0*CONST_PI*CONST_G/(SPEED_C*SPEED_C) * m_densFactor *Dd;

	auto sq = [](double t2, double a_s) {
		return std::sqrt(a_s*a_s + t2);
	};

	auto f = [](double a_s, double sq) {
		return 1.0/(a_s + sq);
	};

	double t2 = theta.getLengthSquared();
	double tx = theta.getX();
	double ty = theta.getY();

	double sq_s = sq(t2, m_scaleRadius);
	double sq_a = sq(t2, m_coreRadius);
	double f_s = f(m_scaleRadius, sq_s);
	double f_a = f(m_coreRadius, sq_a);

	double firstpart = f_s*f_s/sq_s - f_a*f_a/sq_a;
	double secondpart = f_a - f_s;
	axx = Q*( tx*tx*firstpart + secondpart);
	ayy = Q*( ty*ty*firstpart + secondpart);
	axy = Q*( tx*ty*firstpart );
*/

double PIMDLens::getMassInside(double thetaLength) const
{
	double t2 = thetaLength*thetaLength;

	return m_massFactor*(SQRT(m_a2 + t2) - m_coreRadius - SQRT(m_s2 + t2) + m_scaleRadius);
}

double PIMDLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	double t2 = thetaLength*thetaLength;

	return m_densFactor*(1.0/SQRT(m_a2 + t2) - 1.0/SQRT(m_s2 + t2));
}

} // end namespace
