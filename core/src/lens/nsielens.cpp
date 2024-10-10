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
#include "nsielens.h"
#include "constants.h"
#include <iostream>

namespace grale
{

std::unique_ptr<GravitationalLensParams> NSIELensParams::createCopy() const
{
	return std::make_unique<NSIELensParams>(m_velocityDispersion, m_ellipticity, m_angularCoreRadius);
}

bool NSIELensParams::write(serut::SerializationInterface &si) const
{ 
	double array[3]; 
	
	array[0] = m_velocityDispersion; 
	array[1] = m_ellipticity; 
	array[2] = m_angularCoreRadius; 
	
	if (!si.writeDoubles(array, 3)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool NSIELensParams::read(serut::SerializationInterface &si)
{ 
	double array[3]; 
	
	if (!si.readDoubles(array, 3)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	m_velocityDispersion = array[0]; 
	m_ellipticity = array[1]; 
	m_angularCoreRadius = array[2]; 
	
	return true; 
}

NSIELens::NSIELens() : GravitationalLens(GravitationalLens::NSIE)
{
}

NSIELens::~NSIELens()
{
}

bool NSIELens::processParameters(const GravitationalLensParams *params)
{
	const NSIELensParams *p = dynamic_cast<const NSIELensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'NSIELensParams'");
		return false;
	}
	
	V = p->getVelocityDispersion();
	F = p->getEllipticity();
	corerad = p->getAngularCoreRadius();
	
	return true;
}

bool NSIELens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *alpha0) const
{
	double factor = 4.0*CONST_PI*V*V*SQRT(F)/(SPEED_C*SPEED_C*SQRT(1.0-F*F));
	double xfrac = (theta.getX()*SQRT(1.0-F*F))/(SQRT(F*F*theta.getY()*theta.getY()+corerad*corerad+theta.getX()*theta.getX())+F*corerad);
	double yfrac = (F*theta.getY()*SQRT(1.0-F*F))/(F*SQRT(F*F*theta.getY()*theta.getY()+corerad*corerad+theta.getX()*theta.getX())+corerad);

	Vector2D<double> alpha = Vector2D<double>(factor*ATANH(xfrac), factor*ATAN(yfrac));

	*alpha0 = alpha;
	
	return true;
}


bool NSIELens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double factor = 4.0*CONST_PI*V*V*SQRT(F)/(SPEED_C*SPEED_C*SQRT(1.0-F*F));
	double x = theta.getX();
	double y = theta.getY();
	double g = SQRT(F*F*y*y+corerad*corerad+x*x);
	double denom = (F*F*y*y+corerad*corerad+2.0*F*corerad*g+F*F*corerad*corerad+x*x*F*F)*g;
	
	axx = factor*SQRT(1.0-F*F)*(F*F*y*y+corerad*corerad+F*corerad*g)/denom;
	ayy = factor*F*SQRT(1.0-F*F)*(F*corerad*corerad+F*x*x+corerad*g)/denom;
	axy = -factor*x*y*F*F*SQRT(1.0-F*F)/denom;
	
	return true;
}


double NSIELens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double dens;
		
	dens = (V*V/(2.0L*CONST_G*getLensDistance()))*SQRT(F/(theta.getX()*theta.getX()+F*F*theta.getY()*theta.getY()+corerad*corerad));
	
	return dens;
}

bool NSIELens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	Vector2D<double> alpha;

	getAlphaVector(theta, &alpha);

	double factor = 2.0*CONST_PI*V*V*corerad/(SPEED_C*SPEED_C*SQRT(F));
	double part = SQRT(F*F*theta.getY()*theta.getY()+corerad*corerad+theta.getX()*theta.getX())+corerad/F;

	double potential = theta.getX()*alpha.getX()+theta.getY()*alpha.getY() 
		         - factor*LN(part*part+(1.0-F*F)*theta.getY()*theta.getY());
	
	*pPotentialValue = (D_ds/D_s)*potential;
	return true;
}

bool NSIELens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	double factor = 4.0*CONST_PI*V*V*SQRT(F)/(SPEED_C*SPEED_C*SQRT(1.0-F*F));
	*pDeflectionScale = factor;
	*pPotentialScale = factor*factor; // TODO: is this ok?
	return true;
}

bool NSIELens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	if (getenv("GRALE_EXPERIMENTAL_OPENCLPARAMS"))
	{
		*pNumIntParams = 0;
		*pNumFloatParams = 5;
		return true;
	}

	*pNumIntParams = 0;
	*pNumFloatParams = 3;
	return true;
}

bool NSIELens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	if (getenv("GRALE_EXPERIMENTAL_OPENCLPARAMS"))
	{
		std::cerr << "GRALE_EXPERIMENTAL_OPENCLPARAMS on" << std::endl;
		double velScale = 100000.0; // units of 100km/s
		double factor = 4.0*CONST_PI*velScale*velScale/(SPEED_C*SPEED_C);

		pFloatParams[0] = (float)F; // ellipticity
		pFloatParams[1] = (float)(V/velScale); // scaled velocity dispersion
		pFloatParams[2] = (float)(corerad/deflectionScale); // scaled core radius
		pFloatParams[3] = (float)(factor/deflectionScale);
		pFloatParams[4] = (float)(deflectionScale*deflectionScale/potentialScale);
		return true;
	}

	double factor = 4.0*CONST_PI*V*V*SQRT(F)/(SPEED_C*SPEED_C*SQRT(1.0-F*F));
	pFloatParams[0] = (float)F;
	pFloatParams[1] = (float)(corerad/deflectionScale);
	pFloatParams[2] = (float)(factor/deflectionScale);

	return true;
}

std::string NSIELens::getCLProgram(std::string &subRoutineName, bool derivatives, bool potential) const
{
	subRoutineName = "clNSIELensProgram";

	if (getenv("GRALE_EXPERIMENTAL_OPENCLPARAMS"))
	{
		std::string program = R"XYZ(

LensQuantities clNSIELensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
	float F = pFloatParams[0];
	float scaledV = pFloatParams[1];
	float scaledCore = pFloatParams[2];
	float extraFactor = pFloatParams[3];
	float Q = sqrt(1.0f-F*F);

	float denomPart = sqrt(F*F*coord.y*coord.y + scaledCore*scaledCore + coord.x*coord.x);
	float X = (coord.x*Q)/(denomPart + F*scaledCore);
	float Y = (F*coord.y*Q)/(F*denomPart + scaledCore);

	float factor = extraFactor*scaledV*scaledV*sqrt(F)/Q;

	LensQuantities r;
	r.alphaX = factor*atanh(X);
	r.alphaY = factor*atan(Y);
)XYZ";
		if (potential)
			program += R"XYZ(
	float potentialPrefactor = pFloatParams[4];
	float factor2 = 0.5*extraFactor*scaledV*scaledV*scaledCore/sqrt(F);
	float part = sqrt(F*F*coord.y*coord.y + scaledCore*scaledCore + coord.x*coord.x) + scaledCore/F;
	r.potential = (coord.x*r.alphaX + coord.y*r.alphaY
	                       -factor2*log(part*part + (1.0f-F*F)*coord.y*coord.y)
							)*potentialPrefactor;
)XYZ";

		if (derivatives)
			program += R"XYZ(
	float g = sqrt(F*F*coord.y*coord.y + scaledCore*scaledCore + coord.x*coord.x);
	float denom = (F*F*coord.y*coord.y + scaledCore*scaledCore + 2.0f*F*scaledCore*g + F*F*scaledCore*scaledCore + coord.x*coord.x*F*F)*g;
	r.axx = factor*Q*(F*F*coord.y*coord.y + scaledCore*scaledCore + F*scaledCore*g)/denom;
	r.ayy = factor*F*Q*(F*scaledCore*scaledCore + F*coord.x*coord.x + scaledCore*g)/denom;
	r.axy = -factor*coord.x*coord.y*F*F*Q/denom;
)XYZ";

		program += R"XYZ(
	return r;
}
)XYZ";

		return program;
	}

	std::string program = R"XYZ(

LensQuantities clNSIELensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
	float F = pFloatParams[0];
	float w = pFloatParams[1]; // core width
	float factor = pFloatParams[2]; // factor that contains velocity dispersion
	float Q = sqrt(1.0f-F*F);

	float denomPart = sqrt(F*F*coord.y*coord.y + w*w + coord.x*coord.x);
	float X = (coord.x*Q)/(denomPart + F*w);
	float Y = (F*coord.y*Q)/(F*denomPart + w);

	LensQuantities r;
	r.alphaX = factor*atanh(X);
	r.alphaY = factor*atan(Y);

)XYZ";

	if (potential)
		program += R"XYZ(
	r.potential = nan((uint)0); // TODO?
)XYZ";

	if (derivatives)
		program += R"XYZ(
	r.axx = nan((uint)0); // TODO?
	r.ayy = nan((uint)0); // TODO?
	r.axy = nan((uint)0); // TODO?
)XYZ";

	program += R"XYZ(
	return r;
}
	)XYZ";

	return program;
}

} // end namespace

