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

#include "piemdlens.h"
#include "constants.h"
#include "mathfunctions.h"
#include <iostream>

using namespace std;

namespace grale
{

PIEMDLensParams::PIEMDLensParams() 
	: m_sigma0(0),
	  m_coreRadius(0),
	  m_scaleRadius(0),
	  m_epsilon(0)
{
}

PIEMDLensParams::PIEMDLensParams(double sigma0, double coreRadius, double scaleRadius, double epsilon)
	: m_sigma0(sigma0),
	  m_coreRadius(coreRadius),
	  m_scaleRadius(scaleRadius),
	  m_epsilon(epsilon)
{
}

PIEMDLensParams::~PIEMDLensParams()
{
}

bool PIEMDLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_sigma0) &&
	      si.writeDouble(m_coreRadius) &&
		  si.writeDouble(m_scaleRadius) &&
		  si.writeDouble(m_epsilon)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool PIEMDLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_sigma0) &&
		  si.readDouble(&m_coreRadius) &&
		  si.readDouble(&m_scaleRadius) &&
		  si.readDouble(&m_epsilon)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> PIEMDLensParams::createCopy() const
{
	return std::make_unique<PIEMDLensParams>(m_sigma0, m_coreRadius, m_scaleRadius, m_epsilon);
}

PIEMDLens::PIEMDLens() : GravitationalLens(PIEMD)
{
}

PIEMDLens::~PIEMDLens()
{
}

bool PIEMDLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const PIEMDLensParams *pParams = dynamic_cast<const PIEMDLensParams* >(pLensParams);
	if (!pParams)
	{
		setErrorString("Parameters are not of type 'PIEMDLensParams'");
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

	m_epsilon = pParams->getEpsilon();
	if (m_epsilon < 0.0 || m_epsilon >= 1.0)
	{
		setErrorString("Invalid epsilon value");
		return false;
	}

	m_sqrtEpsilon = SQRT(m_epsilon);
	m_epsFrac = (1.0-m_epsilon)/(1.0+m_epsilon);

	double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());
	
	m_deflectionFactor = m_sigma0/densScale * (1.0-m_epsilon*m_epsilon)/SQRT(m_epsilon);
	m_angularScale = m_scaleRadius/2.0;

	m_a2 = m_coreRadius*m_coreRadius;
	m_s2 = m_scaleRadius*m_scaleRadius;
	m_densFactor = m_sigma0*m_scaleRadius*m_coreRadius/(m_scaleRadius-m_coreRadius);

	return true;
}

Vector2Dd PIEMDLens::calcI(double omega, Vector2Dd theta) const
{
	double X = theta.getX()/m_angularScale;
	double Y = theta.getY()/m_angularScale;
	double W = omega/m_angularScale;

	double A = X*m_epsFrac;
	double tmpX = X/(1.0+m_epsilon);
	double tmpY = Y/(1.0-m_epsilon);
	double B = -Y/m_epsFrac + 2.0*m_sqrtEpsilon*SQRT(W*W + tmpX*tmpX + tmpY*tmpY);
	double C = X;
	double D = -Y+2.0*W*m_sqrtEpsilon;
	double C2D2 = C*C+D*D;

	double H = (B*C-D*A)/C2D2;
	double G = (C*A+D*B)/C2D2;

	double re = atan2(H, G) * omega;
	double im = -LN(SQRT(G*G+H*H)) * omega;
	return Vector2Dd(re, im);
}

bool PIEMDLens::getAlphaVector(Vector2Dd theta,Vector2Dd *pAlpha) const
{
	Vector2Dd corePart = calcI(m_coreRadius, theta);
	Vector2Dd scalePart = calcI(m_scaleRadius, theta);

	double diffRadii = m_scaleRadius-m_coreRadius;
	*pAlpha = (corePart*(m_scaleRadius/diffRadii) - scalePart*(m_coreRadius/diffRadii))*m_deflectionFactor;
	return true;
}

void PIEMDLens::getIDerivs(double omega, Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double X = theta.getX()/m_angularScale;
	double Y = theta.getY()/m_angularScale;
	double W = omega/m_angularScale;
	double tmpX = X/(1.0+m_epsilon);
	double tmpY = Y/(1.0-m_epsilon);

	double I0 = m_deflectionFactor;
	double Q = m_epsFrac * X;
	double sqrtOmegaEtc = SQRT(W*W+tmpX*tmpX+tmpY*tmpY);
	double P = Y/m_epsFrac -2.0*m_sqrtEpsilon*sqrtOmegaEtc;
	double U = X;
	double V = Y-2.0*W*m_sqrtEpsilon;

	double U2V2 = U*U+V*V;
	double R = (Q*U+P*V)/U2V2;
	double Z = (Q*V-P*U)/U2V2;

	double Qx = m_epsFrac;
	double Px = -2.0*m_sqrtEpsilon/sqrtOmegaEtc * tmpX/(1.0+m_epsilon);
	double Py = 1.0/m_epsFrac - 2.0*m_sqrtEpsilon/sqrtOmegaEtc * tmpY/(1.0-m_epsilon);

	double U2V2Squared = U2V2*U2V2;
	double QUPVoverU2V22 = (Q*U+P*V)/U2V2Squared;
	double QVPUoverU2V22 = (Q*V-P*U)/U2V2Squared;
	double Rx = (Qx*U + Q + Px*V)/U2V2 - QUPVoverU2V22 * 2.0*U;
	double Ry = (Py*V + P)/U2V2 - QUPVoverU2V22 * 2.0*V;
	double Zx = (Qx*V - Px*U - P)/U2V2 - QVPUoverU2V22 * 2.0*U;
	double Zy = (Q - Py*U)/U2V2 - QVPUoverU2V22 * 2.0*V;

	double R2Z2 = R*R+Z*Z;
	double factor = I0*omega/R2Z2*1.0/m_angularScale;

	axx = factor*(R*Zx - Z*Rx);
	ayy = -factor*(R*Ry + Z*Zy);
	axy = -factor*(R*Rx + Z*Zx);
}

bool PIEMDLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	//return GravitationalLens::getAlphaVectorDerivatives(theta, axx, ayy, axy);

	double axx_a, ayy_a, axy_a;
	double axx_s, ayy_s, axy_s;

	getIDerivs(m_coreRadius, theta, axx_a, ayy_a, axy_a);
	getIDerivs(m_scaleRadius, theta, axx_s, ayy_s, axy_s);

	double diffRadii = m_scaleRadius-m_coreRadius;
	double fa = m_scaleRadius/diffRadii;
	double fs = m_coreRadius/diffRadii;
	
	axx = axx_a*fa - axx_s*fs;
	ayy = ayy_a*fa - ayy_s*fs;
	axy = axy_a*fa - axy_s*fs;
	return true;
}

double PIEMDLens::getSurfaceMassDensity(Vector2Dd theta) const
{
	double tmpX = theta.getX()/(1.0+m_epsilon);
	double tmpY = theta.getY()/(1.0-m_epsilon);

	double thetaLength = SQRT(tmpX*tmpX + tmpY*tmpY);
	double t2 = thetaLength*thetaLength;

	return m_densFactor*(1.0/SQRT(m_a2 + t2) - 1.0/SQRT(m_s2 + t2));
}

} // end namespace

