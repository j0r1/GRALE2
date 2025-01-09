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
#include "utils.h"
#include <iostream>

using namespace std;

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

	m_angularScale = 8.0*CONST_PI*CONST_G/(SPEED_C*SPEED_C) * m_densFactor * Dd;

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

bool PIMDLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	*pDeflectionScale = m_angularScale;
	*pPotentialScale = m_angularScale*m_angularScale;
	return true;
}

bool PIMDLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 3;
	return true;
}

bool PIMDLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());

	pFloatParams[0] = (float)(m_sigma0/densScale);
	pFloatParams[1] = (float)(m_coreRadius/deflectionScale);
	pFloatParams[2] = (float)(m_scaleRadius/deflectionScale);
	return true;
}

string PIMDLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	subRoutineName = "clPIMDLensProgram";

	double Q_scaled = m_angularScale/deflectionScale;

	string pimdParams = R"XYZ(
	float sigma0 = pFloatParams[0];
	float coreRadius = pFloatParams[1];
	float scaleRadius = pFloatParams[2];
	)XYZ";

	string program = R"XYZ(

LensQuantities clPIMDLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{)XYZ" + pimdParams + R"XYZ(

	float Q_scaled = )XYZ" + float_to_string((float)Q_scaled) + R"XYZ(;
	float t2 = coord.x*coord.x + coord.y*coord.y;

	double sq_s = sqrt(scaleRadius*scaleRadius + t2);
	double sq_a = sqrt(coreRadius*coreRadius + t2);
	double f_s = 1.0/(scaleRadius + sq_s);
	double f_a = 1.0/(coreRadius + sq_a);
	float f_diff = f_a - f_s;
	float common = Q_scaled*f_diff;

	LensQuantities r;
	r.alphaX = common*coord.x;
	r.alphaY = common*coord.y;
)XYZ";
	if (potential)
		program += R"XYZ(
	r.potential = nan((uint)0); // Currently not supported
)XYZ";
	if (derivatives)
		program += R"XYZ(
	float firstpart = f_s*f_s/sq_s - f_a*f_a/sq_a;
	float secondpart = f_diff;
	float tx = coord.x;
	float ty = coord.y;
	r.axx = Q_scaled*( tx*tx*firstpart + secondpart);
	r.ayy = Q_scaled*( ty*ty*firstpart + secondpart);
	r.axy = Q_scaled*( tx*ty*firstpart );
)XYZ";

	program += R"XYZ(
	return r;

}
	)XYZ";

	return program;
}

unique_ptr<GravitationalLensParams> PIMDLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());

	double sigma0 = (double)pFloatParams[0] * densScale;
	double coreRad = (double)pFloatParams[1] * deflectionScale;
	double scaleRad = (double)pFloatParams[2] * deflectionScale;

	return make_unique<PIMDLensParams>(sigma0, coreRad, scaleRad);
}

vector<CLFloatParamInfo> PIMDLens::getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const
{
	double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());

	return {
		{ .name = "centraldensity_scaled", .offset = 0, .scaleFactor = densScale, .hardMin = 0 },
		{ .name = "coreradius_scaled", .offset = 1, .scaleFactor = deflectionScale, .hardMin = 0 },
		{ .name = "scaleradius_scaled", .offset = 2, .scaleFactor = deflectionScale, .hardMin = 0 },
	};
}

} // end namespace
