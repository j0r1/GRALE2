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
#include "utils.h"
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

PIEMDLens::PIEMDLens(LensType t) : GravitationalLens(t)
{
}

PIEMDLens::~PIEMDLens()
{
}

bool PIEMDLens::subInit(double sigma0, double coreRad, double scaleRad, double eps)
{
	m_sigma0 = sigma0;
	//pParams->getCentralDensity();

	m_coreRadius = coreRad; //pParams->getCoreRadius();
	m_scaleRadius = scaleRad; //pParams->getScaleRadius();
	if (m_coreRadius < 0)
	{
		setErrorString("Invalid core radius");
		return false;
	}

	if (m_scaleRadius <= m_coreRadius)
		std::cerr << "WARNING: scale radius <= core radius, lens effect will be disabled" << endl;

	m_epsilon = eps; //pParams->getEpsilon();
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

bool PIEMDLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const PIEMDLensParams *pParams = dynamic_cast<const PIEMDLensParams* >(pLensParams);
	if (!pParams)
	{
		setErrorString("Parameters are not of type 'PIEMDLensParams'");
		return false;
	}

	return subInit(pParams->getCentralDensity(), pParams->getCoreRadius(), pParams->getScaleRadius(), pParams->getEpsilon());
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
	if (m_scaleRadius <= m_coreRadius)
	{
		*pAlpha = Vector2Dd(0,0);
		return true;
	}

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
	if (m_scaleRadius <= m_coreRadius)
	{
		axx = 0;
		ayy = 0;
		axy = 0;
		return true;
	}

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
	if (m_scaleRadius <= m_coreRadius)
		return 0;

	double tmpX = theta.getX()/(1.0+m_epsilon);
	double tmpY = theta.getY()/(1.0-m_epsilon);

	double thetaLength = SQRT(tmpX*tmpX + tmpY*tmpY);
	double t2 = thetaLength*thetaLength;

	return m_densFactor*(1.0/SQRT(m_a2 + t2) - 1.0/SQRT(m_s2 + t2));
}

bool PIEMDLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	*pDeflectionScale = m_angularScale;
	*pPotentialScale = m_angularScale*m_angularScale;
	return true;
}

bool PIEMDLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 4;
	return true;
}

bool PIEMDLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());

	pFloatParams[0] = (float)(m_sigma0/densScale);
	pFloatParams[1] = (float)(m_coreRadius/deflectionScale);
	pFloatParams[2] = (float)(m_scaleRadius/deflectionScale);
	pFloatParams[3] = (float)m_epsilon;
	return true;
}

static inline string getPIEMDSubRoutines(const string &subRoutineName)
{
	string program = R"XYZ(

float2 )XYZ" + subRoutineName + R"XYZ(_calcI(float omega_scaled, float2 theta, float epsilon)
{
	float X = theta.x;
	float Y = theta.y;
	float W = omega_scaled;
	float epsFrac = (1.0f-epsilon)/(1.0f+epsilon);

	float A = X*epsFrac;
	float tmpX = X/(1.0f+epsilon);
	float tmpY = Y/(1.0f-epsilon);
	float sqrtEpsilon = sqrt(epsilon);
	float B = -Y/epsFrac + 2.0f*sqrtEpsilon*sqrt(W*W + tmpX*tmpX + tmpY*tmpY);
	float C = X;
	float D = -Y+2.0f*W*sqrtEpsilon;
	float C2D2 = C*C+D*D;

	float H = (B*C-D*A)/C2D2;
	float G = (C*A+D*B)/C2D2;

	float re = atan2(H, G) * omega_scaled;
	float im = -log(sqrt(G*G+H*H)) * omega_scaled;
	return (float2)(re, im);
}

void )XYZ" + subRoutineName + R"XYZ(_calcIDerivs(float omega_scaled, float2 theta, float *pAxx, float *pAyy, float *pAxy, float epsilon, float deflectionFactor)
{
	float X = theta.x;
	float Y = theta.y;
	float W = omega_scaled;
	float epsFrac = (1.0f-epsilon)/(1.0f+epsilon);
	float tmpX = X/(1.0f+epsilon);
	float tmpY = Y/(1.0f-epsilon);
	float sqrtEpsilon = sqrt(epsilon);

	float I0 = deflectionFactor;
	float Q = epsFrac * X;
	float sqrtOmegaEtc = sqrt(W*W+tmpX*tmpX+tmpY*tmpY);
	float P = Y/epsFrac -2.0f*sqrtEpsilon*sqrtOmegaEtc;
	float U = X;
	float V = Y-2.0f*W*sqrtEpsilon;

	float U2V2 = U*U+V*V;
	float R = (Q*U+P*V)/U2V2;
	float Z = (Q*V-P*U)/U2V2;

	float Qx = epsFrac;
	float Px = -2.0f*sqrtEpsilon/sqrtOmegaEtc * tmpX/(1.0f+epsilon);
	float Py = 1.0/epsFrac - 2.0f*sqrtEpsilon/sqrtOmegaEtc * tmpY/(1.0f-epsilon);

	float U2V2Squared = U2V2*U2V2;
	float QUPVoverU2V22 = (Q*U+P*V)/U2V2Squared;
	float QVPUoverU2V22 = (Q*V-P*U)/U2V2Squared;
	float Rx = (Qx*U + Q + Px*V)/U2V2 - QUPVoverU2V22 * 2.0f*U;
	float Ry = (Py*V + P)/U2V2 - QUPVoverU2V22 * 2.0f*V;
	float Zx = (Qx*V - Px*U - P)/U2V2 - QVPUoverU2V22 * 2.0f*U;
	float Zy = (Q - Py*U)/U2V2 - QVPUoverU2V22 * 2.0f*V;

	float R2Z2 = R*R+Z*Z;
	float factor = I0*omega_scaled/R2Z2;

	*pAxx = factor*(R*Zx - Z*Rx);
	*pAyy = -factor*(R*Ry + Z*Zy);
	*pAxy = -factor*(R*Rx + Z*Zx);
}
)XYZ";
	return program;
}

static inline string getPIEMDMainProgram(const string &subRoutineName, const string &piemdParams, bool derivatives, bool potential)
{
	string program = getPIEMDSubRoutines(subRoutineName);
	program += R"XYZ(

LensQuantities )XYZ" + subRoutineName + R"XYZ((float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
)XYZ";

	program += piemdParams;

	program += R"XYZ(
	LensQuantities r;

	float2 corePart = )XYZ" + subRoutineName + R"XYZ(_calcI(coreRadius, coord, epsilon);
	float2 scalePart = )XYZ" + subRoutineName + R"XYZ(_calcI(scaleRadius, coord, epsilon);
	float deflectionFactor = sigma0 * (1.0f-epsilon*epsilon)/sqrt(epsilon);

	float diffRadii = scaleRadius-coreRadius;
	float fa = 0;
	float fs = 0;
	if (diffRadii <= 0)
	{
		r.alphaX = 0;
		r.alphaY = 0;
	}
	else
	{
		fa = scaleRadius/diffRadii;
		fs = coreRadius/diffRadii;

		float2 alpha = (corePart*fa - scalePart*fs) * deflectionFactor;

		//if (isnan(alpha.x) || isnan(alpha.y))
		//	printf("NaN in alpha for point %d, genome %d\n", get_global_id(0), get_global_id(1));

		r.alphaX = alpha.x;
		r.alphaY = alpha.y;
	}
)XYZ";
	if (potential)
		program += R"XYZ(
	r.potential = nan((uint)0); // Currently not supported
)XYZ";
	if (derivatives)
		program += R"XYZ(
	if (diffRadii <= 0) // this means coreRadius == scaleRadius, and a.._a == a.._s
	{
		r.axx = 0;
		r.ayy = 0;
		r.axy = 0;
	}
	else
	{
		float axx_a, ayy_a, axy_a;
		float axx_s, ayy_s, axy_s;

		)XYZ" + subRoutineName + R"XYZ(_calcIDerivs(coreRadius, coord, &axx_a, &ayy_a, &axy_a, epsilon, deflectionFactor);
		)XYZ" + subRoutineName + R"XYZ(_calcIDerivs(scaleRadius, coord, &axx_s, &ayy_s, &axy_s, epsilon, deflectionFactor);

		r.axx = axx_a*fa - axx_s*fs;
		r.ayy = ayy_a*fa - ayy_s*fs;
		r.axy = axy_a*fa - axy_s*fs;
	}

	//if (isnan(r.axx) || isnan(r.ayy) || isnan(r.axy))
	//	printf("NaN in alpha derivs for point %d, genome %d\n", get_global_id(0), get_global_id(1));

)XYZ";
	program += R"XYZ(
	return r;
}
	)XYZ";

	return program;
}

string PIEMDLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	subRoutineName = "clPIEMDLensProgram";

	string piemdParams = R"XYZ(
	float sigma0 = pFloatParams[0];
	float coreRadius = pFloatParams[1];
	float scaleRadius = pFloatParams[2];
	float epsilon = pFloatParams[3];
	)XYZ";

	return getPIEMDMainProgram(subRoutineName, piemdParams, derivatives, potential);
}

unique_ptr<GravitationalLensParams> PIEMDLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());

	double sigma0 = (double)pFloatParams[0] * densScale;
	double coreRad = (double)pFloatParams[1] * deflectionScale;
	double scaleRad = (double)pFloatParams[2] * deflectionScale;
	double eps = (double)pFloatParams[3];

	return make_unique<PIEMDLensParams>(sigma0, coreRad, scaleRad, eps);
}

vector<CLFloatParamInfo> PIEMDLens::getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const
{
	double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());

	return {
		{ .name = "centraldensity_scaled", .offset = 0, .scaleFactor = densScale, .hardMin = 0 },
		{ .name = "coreradius_scaled", .offset = 1, .scaleFactor = deflectionScale, .hardMin = 0 },
		{ .name = "scaleradius_scaled", .offset = 2, .scaleFactor = deflectionScale, .hardMin = 0 },
		{ .name = "epsilon", .offset = 3, .scaleFactor = 1.0, .hardMin = 0.01, .hardMax = 0.99 }, // TODO 0 and 1 are not allowed, what are good bounds?
	};
}

LTPIEMDLensParams::LTPIEMDLensParams() 
	: m_velDisp(0),
	  m_coreRadius(0),
	  m_scaleRadius(0),
	  m_ellipticity(0)
{
}

LTPIEMDLensParams::LTPIEMDLensParams(double velDisp, double coreRadius, double scaleRadius, double ellipticity)
	: m_velDisp(velDisp),
	  m_coreRadius(coreRadius),
	  m_scaleRadius(scaleRadius),
	  m_ellipticity(ellipticity)
{
}

LTPIEMDLensParams::~LTPIEMDLensParams()
{
}

bool LTPIEMDLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_velDisp) &&
	      si.writeDouble(m_coreRadius) &&
		  si.writeDouble(m_scaleRadius) &&
		  si.writeDouble(m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool LTPIEMDLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_velDisp) &&
		  si.readDouble(&m_coreRadius) &&
		  si.readDouble(&m_scaleRadius) &&
		  si.readDouble(&m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> LTPIEMDLensParams::createCopy() const
{
	return std::make_unique<LTPIEMDLensParams>(m_velDisp, m_coreRadius, m_scaleRadius, m_ellipticity);
}

LTPIEMDLens::LTPIEMDLens() : PIEMDLens(LTPIEMD)
{
}

LTPIEMDLens::~LTPIEMDLens()
{
}

bool LTPIEMDLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const LTPIEMDLensParams *pParams = dynamic_cast<const LTPIEMDLensParams* >(pLensParams);
	if (!pParams)
	{
		setErrorString("Parameters are not of type 'LTPIEMDLensParams'");
		return false;
	}

	double velDisp = pParams->getVelocityDispersion();
	double epsHat = pParams->getEllipticity();
	double a = pParams->getCoreRadius();
	double s = pParams->getScaleRadius();
	double Dd = getLensDistance();

	double centralDens = (3.0*velDisp*velDisp)/(4.0*CONST_G*Dd)*(s*s-a*a)/(a*s*s);
	double e = 1.0-std::sqrt((1.0-epsHat)/(1.0+epsHat));
	double eps = e/(2.0-e);

	m_velDisp = velDisp;
	m_epsHat = epsHat;

	return subInit(centralDens, a, s, eps);
}

static const double velScale = 100000.0; // units of 100km/s

bool LTPIEMDLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	pFloatParams[0] = (float)(m_velDisp/velScale);
	pFloatParams[1] = (float)(m_coreRadius/deflectionScale);
	pFloatParams[2] = (float)(m_scaleRadius/deflectionScale);
	pFloatParams[3] = (float)m_epsHat;
	return true;
}

string LTPIEMDLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	subRoutineName = "clLTPIEMDLensProgram";

	//double densScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance()); // sigma0 should be in these units
	//double sigma0ScaleFactor = (velScale*velScale)/(4.0*CONST_G*Dd*deflectionScale);
	//double finalDensScale = sigma0ScaleFactor/densScale;
	double finalDensScale = ((velScale/SPEED_C)*(velScale/SPEED_C)*CONST_PI)/deflectionScale;

	string piemdParams = R"XYZ(
	float velDisp = pFloatParams[0]; // scaled with 'velScale'
	float coreRadius = pFloatParams[1];
	float scaleRadius = pFloatParams[2];
	float epsHat = pFloatParams[3];
	
	float finalDensScale = )XYZ" + float_to_string((float)finalDensScale) + R"XYZ(;
	float aRescaled = coreRadius/scaleRadius;
	float sigma0 = finalDensScale * (3.0*velDisp*velDisp) * (1.0 - aRescaled*aRescaled)/coreRadius;
	float eTmp = 1.0-sqrt((1.0-epsHat)/(1.0+epsHat));
	float epsilon = eTmp/(2.0-eTmp);
	)XYZ";

	return getPIEMDMainProgram(subRoutineName, piemdParams, derivatives, potential);
}

unique_ptr<GravitationalLensParams> LTPIEMDLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	double velDisp = (double)pFloatParams[0] * velScale;
	double coreRad = (double)pFloatParams[1] * deflectionScale;
	double scaleRad = (double)pFloatParams[2] * deflectionScale;
	double epsHat = (double)pFloatParams[3];

	return make_unique<LTPIEMDLensParams>(velDisp, coreRad, scaleRad, epsHat);
}

vector<CLFloatParamInfo> LTPIEMDLens::getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const
{
	return {
		{ .name = "velocitydispersion_scaled", .offset = 0, .scaleFactor = velScale, .hardMin = 0 },
		{ .name = "coreradius_scaled", .offset = 1, .scaleFactor = deflectionScale, .hardMin = 0 },
		{ .name = "scaleradius_scaled", .offset = 2, .scaleFactor = deflectionScale, .hardMin = 0 },
		{ .name = "ellipticity", .offset = 3, .scaleFactor = 1.0, .hardMin = 0.01, .hardMax = 0.99 }, // TODO 0 and 1 are not allowed, what are good bounds?
	};
}

} // end namespace

