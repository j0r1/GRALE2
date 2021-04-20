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

#include "alphapotlens.h"
#include "constants.h"
#include "mathfunctions.h"
#include <iostream>

using namespace std;

namespace grale
{

AlphaPotLensParams::AlphaPotLensParams() 
	: m_b(0),
	  m_s(0),
	  m_q(1),
	  m_K2(0),
	  m_alpha(0)
{
}

AlphaPotLensParams::AlphaPotLensParams(double b, double s, double q, double K2, double alpha)
	: m_b(b),
	  m_s(s),
	  m_q(q),
	  m_K2(K2),
	  m_alpha(alpha)
{
}

AlphaPotLensParams::~AlphaPotLensParams()
{
}

bool AlphaPotLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_b) &&
	      si.writeDouble(m_s) &&
		  si.writeDouble(m_q) &&
		  si.writeDouble(m_K2) &&
		  si.writeDouble(m_alpha)
		  ))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool AlphaPotLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_b) &&
		  si.readDouble(&m_s) &&
		  si.readDouble(&m_q) &&
		  si.readDouble(&m_K2) &&
		  si.readDouble(&m_alpha)
		  ))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> AlphaPotLensParams::createCopy() const
{
	return std::make_unique<AlphaPotLensParams>(m_b, m_s, m_q, m_K2, m_alpha);
}

AlphaPotLens::AlphaPotLens() : GravitationalLens(AlphaPot)
{
}

AlphaPotLens::~AlphaPotLens()
{
}

bool AlphaPotLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const AlphaPotLensParams *pParams = dynamic_cast<const AlphaPotLensParams* >(pLensParams);
	if (!pParams)
	{
		setErrorString("Parameters are not of type 'AlphaPotLensParams'");
		return false;
	}

	m_b = pParams->getB();
	m_s2 = pParams->getS() * pParams->getS();
	m_q2 = pParams->getQ() * pParams->getQ();
	m_K2 = pParams->getK2();
	m_alpha = pParams->getAlpha();

	if (m_q2 == 0)
	{
		setErrorString("Invalid q value");
		return false;
	}

	return true;
}

bool AlphaPotLens::getAlphaVector(Vector2Dd theta,Vector2Dd *pAlpha) const
{
	double x = theta.getX();
	double y = theta.getY();
	double factor = 0.5*m_b*m_alpha*POW(m_s2 + x*x + y*y/m_q2 + m_K2*x*y, m_alpha*0.5-1.0);
	double ax = (2.0*x+m_K2*y)*factor;
	double ay = (2.0*y/m_q2 + m_K2*x)*factor;
	*pAlpha = Vector2Dd(ax, ay);
	return true;
}

bool AlphaPotLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	// return GravitationalLens::getAlphaVectorDerivatives(theta, axx, ayy, axy);
	double x = theta.getX();
	double y = theta.getY();
	double B = m_s2 + x*x + y*y/m_q2 + m_K2*x*y;
	double A = 0.5*m_b*m_alpha*POW(B, m_alpha*0.5-2.0);
	double xfactor = 2.0*x+m_K2*y;
	double yfactor = 2.0*y/m_q2 + m_K2*x;
	double halfAlphaMinusOne = (0.5*m_alpha-1.0);

	axx = A*(halfAlphaMinusOne*xfactor*xfactor+2.0*B);
	ayy = A*(halfAlphaMinusOne*yfactor*yfactor+(2.0/m_q2)*B);
	axy = A*(halfAlphaMinusOne*xfactor*yfactor+m_K2*B);
	return true;
}

double AlphaPotLens::getSurfaceMassDensity(Vector2Dd theta) const
{
	double x = theta.getX();
	double y = theta.getY();
	double B = m_s2 + x*x + y*y/m_q2 + m_K2*x*y;
	double A = 0.25*m_b*m_alpha*POW(B, m_alpha*0.5-2.0);
	double xfactor = 2.0*x+m_K2*y;
	double yfactor = 2.0*y/m_q2 + m_K2*x;
	double halfAlphaMinusOne = (0.5*m_alpha-1.0);

	double sigmaScale = (SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G))/getLensDistance();
	return sigmaScale*A*(
			halfAlphaMinusOne*(xfactor*xfactor+yfactor*yfactor) +
			(2.0+2.0/m_q2)*B
			);
}

bool AlphaPotLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	double x = theta.getX();
	double y = theta.getY();
	double phi = m_b*POW(m_s2 + x*x + y*y/m_q2 + m_K2*x*y, m_alpha*0.5);

	*pPotentialValue = D_ds/D_s * phi;
	return true;
}

} // end namespace

