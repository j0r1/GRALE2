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
#include "plummerlens.h"
#include "constants.h"
#include "utils.h"
#include <iostream>

namespace grale
{

std::unique_ptr<GravitationalLensParams> PlummerLensParams::createCopy() const
{
	return std::make_unique<PlummerLensParams>(lensmass,angwidth);
}

bool PlummerLensParams::write(serut::SerializationInterface &si) const
{ 
	double array[2]; 
	
	array[0] = lensmass; 
	array[1] = angwidth; 
	
	if (!si.writeDoubles(array,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool PlummerLensParams::read(serut::SerializationInterface &si)
{ 
	double array[2]; 
	
	if (!si.readDoubles(array,2)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	lensmass = array[0]; 
	angwidth = array[1]; 
	
	return true; 
}

PlummerLens::PlummerLens() : SymmetricLens(GravitationalLens::Plummer)
{
}

PlummerLens::~PlummerLens()
{
}

bool PlummerLens::processParameters(const GravitationalLensParams *params)
{
	const PlummerLensParams *p = dynamic_cast<const PlummerLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'PlummerLensParams'");
		return false;
	}
	
	mass = p->getLensMass();
	angularwidth = p->getAngularWidth();
	angularwidth2 = angularwidth*angularwidth;
	return true;
}

double PlummerLens::getMassInside(double thetaLength) const
{
	double t = thetaLength/angularwidth;
	double massinside;
	double x = t*t;
	double x2 = 1.0/x;

	massinside = mass*(1.0/(1.0+x2));
	return massinside;
}

double PlummerLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	double scaledtheta = thetaLength/angularwidth;
	double denom = (1.0+scaledtheta*scaledtheta)*getLensDistance();
	double dens =  mass/(CONST_PI*angularwidth2*denom*denom);
	return dens;
}

bool PlummerLens::getSurfaceMassDensityDerivative(double thetaLength, double &deriv) const
{
	double scaledTheta = thetaLength/angularwidth;
	double denom1 = scaledTheta*scaledTheta + 1.0;
	double denom3 = denom1*denom1*denom1;
	double angularWidth2 = angularwidth*angularwidth;
	double angularWidth3 = angularWidth2*angularwidth;
	double Dd = getLensDistance();
	deriv = (-4.0*mass/(CONST_PI*Dd*Dd*angularWidth3))*(scaledTheta/denom3);
	return true;
}

bool PlummerLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
	                                double *pPotentialValue) const
{
	Vector2D<double> scaledTheta = theta/angularwidth;

	*pPotentialValue = ((D_ds/D_s)*(4.0*CONST_G*mass)/(SPEED_C*SPEED_C*getLensDistance()))*0.5*LN(scaledTheta.getLengthSquared() + 1);
	return true;
}

bool PlummerLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	*pDeflectionScale = SQRT((4.0*CONST_G*mass)/(SPEED_C*SPEED_C*getLensDistance()));
	*pPotentialScale = (4.0*CONST_G*mass)/(SPEED_C*SPEED_C*getLensDistance());
	return true;
}

bool PlummerLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 2;
	return true;
}

bool PlummerLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	// mass factor
	pFloatParams[0] = (float)((4.0*CONST_G*mass)/(SPEED_C*SPEED_C*getLensDistance()*deflectionScale*deflectionScale));
	// width factor
	pFloatParams[1] = (float)(angularwidth/deflectionScale);
	// prefactor needed for potential, to express in potentialscale units;  moved as a constant to the code
	// pFloatParams[2] = (float)(deflectionScale*deflectionScale/(potentialScale*2.0));
	return true;
}

std::unique_ptr<GravitationalLensParams> PlummerLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	float massfac = pFloatParams[0];
	float scaledwidth = pFloatParams[1];
	double mass = ((double)massfac)*SPEED_C*SPEED_C*getLensDistance()*deflectionScale*deflectionScale/(4.0*CONST_G);
	double width = scaledwidth*deflectionScale;

	return std::make_unique<PlummerLensParams>(mass, width);
}

std::vector<CLFloatParamInfo> PlummerLens::getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const
{
	return { 
		{ .name = "mass_scaled", .offset = 0, .scaleFactor = SPEED_C*SPEED_C*getLensDistance()*deflectionScale*deflectionScale/(4.0*CONST_G), .hardMin = 0 },
		{ .name = "width_scaled", .offset = 1, .scaleFactor = deflectionScale, .hardMin = 0 }
	};
}

std::string PlummerLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	std::string program;
	float potentialPrefactor =(float)(deflectionScale*deflectionScale/(potentialScale*2.0));

	program = R"XYZ(
LensQuantities clPlummerLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
	float scaledMass = pFloatParams[0];
	float scaledWidth = pFloatParams[1];
	float scaledWidth2 = scaledWidth*scaledWidth;
	float denom = coord.x*coord.x + coord.y*coord.y + scaledWidth2;
	float denom2 = denom*denom;
	float factor = scaledMass/denom;
	float factor2 = scaledMass/denom2;

	LensQuantities r;
	r.alphaX = factor*coord.x;
	r.alphaY = factor*coord.y;
)XYZ";
	if (potential)
		program += R"XYZ(
	float potentialPrefactor = (float))XYZ" + float_to_string(potentialPrefactor) + R"XYZ(;
	r.potential = potentialPrefactor*scaledMass*log(denom);
)XYZ";
	if (derivatives)
		program += R"XYZ(
	r.axx = factor2*(-coord.x*coord.x+coord.y*coord.y+scaledWidth2);
	r.ayy = factor2*(-coord.y*coord.y+coord.x*coord.x+scaledWidth2);
	r.axy = factor2*(-2.0*coord.x*coord.y);
)XYZ";

	program += R"XYZ(
	return r;
}
)XYZ";
	subRoutineName = "clPlummerLensProgram";
	return program;
}

} // end namespace

