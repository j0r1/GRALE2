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
	double massinside,r2;
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
	*pNumFloatParams = 3;
	return true;
}

bool PlummerLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	double scaledWidth = angularwidth/deflectionScale;

	pFloatParams[0] = (float)((4.0*CONST_G*mass)/(SPEED_C*SPEED_C*getLensDistance()*deflectionScale*deflectionScale));
	pFloatParams[1] = (float)(scaledWidth*scaledWidth);
	pFloatParams[2] = (float)(((2.0*CONST_G*mass)/(SPEED_C*SPEED_C*getLensDistance()))/potentialScale);
	return true;
}

std::string PlummerLens::getCLProgram(std::string &subRoutineName) const
{
	std::string program;

	program += "LensQuantities clPlummerLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)\n";
	program += "{\n";
	program += "	float denom = coord.x*coord.x+coord.y*coord.y+pFloatParams[1];\n";
	program += "	float denom2 = denom*denom;\n";
	program += "	float factor = pFloatParams[0]/denom;\n";
	program += "	float factor2 = pFloatParams[0]/denom2;\n";
	program += "\n";
	program += "	LensQuantities r;\n";
	program += "	r.alphaX = factor*coord.x;\n";
	program += "	r.alphaY = factor*coord.y;\n";
	program += "	r.potential = pFloatParams[2]*log(denom);\n";
	program += "	r.axx = factor2*(-coord.x*coord.x+coord.y*coord.y+pFloatParams[1]);\n";
	program += "	r.ayy = factor2*(-coord.y*coord.y+coord.x*coord.x+pFloatParams[1]);\n";
	program += "	r.axy = factor2*(-2.0*coord.x*coord.y);\n";
	program += "	return r;\n";
	program += "}\n";

	subRoutineName = "clPlummerLensProgram";

	return program;
}

} // end namespace

