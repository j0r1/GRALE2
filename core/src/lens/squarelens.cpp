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
#include "squarelens.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> SquareLensParams::createCopy() const
{
	return std::make_unique<SquareLensParams>(lensmass,angwidth);
}

bool SquareLensParams::write(serut::SerializationInterface &si) const
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

bool SquareLensParams::read(serut::SerializationInterface &si)
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

SquareLens::SquareLens() : GravitationalLens(GravitationalLens::Square)
{
}

SquareLens::~SquareLens()
{
}

// TODO: check if this works for negative masses!!!
bool SquareLens::processParameters(const GravitationalLensParams *params)
{
	const SquareLensParams *p = dynamic_cast<const SquareLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'SquareLensParams'");
		return false;
	}
	
	mass = p->getLensMass();
	angularwidth = p->getAngularWidth();
	angularwidth2 = angularwidth*angularwidth;
	angularscale2 = ((4.0*CONST_G*ABS(mass))/(SPEED_C*SPEED_C*getLensDistance()));
	angularscale = SQRT(angularscale2);
	relativeangularwidth2 = angularwidth2/angularscale2;
	rawd2 = SQRT(relativeangularwidth2)/2.0;
	dens = mass/(getLensDistance()*getLensDistance()*angularwidth2);
	return true;
}

bool SquareLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *alpha) const
{
	Vector2D<double> scaledtheta = theta/angularscale;
	double ax = deflectionFunction1(scaledtheta.getX()-rawd2,scaledtheta.getY()-rawd2)
		  - deflectionFunction1(scaledtheta.getX()+rawd2,scaledtheta.getY()-rawd2)
		  - deflectionFunction1(scaledtheta.getX()-rawd2,scaledtheta.getY()+rawd2)
		  + deflectionFunction1(scaledtheta.getX()+rawd2,scaledtheta.getY()+rawd2);
	double ay = deflectionFunction2(scaledtheta.getX()-rawd2,scaledtheta.getY()-rawd2)
		  - deflectionFunction2(scaledtheta.getX()+rawd2,scaledtheta.getY()-rawd2)
		  - deflectionFunction2(scaledtheta.getX()-rawd2,scaledtheta.getY()+rawd2)
		  + deflectionFunction2(scaledtheta.getX()+rawd2,scaledtheta.getY()+rawd2);

	Vector2D<double> a(ax,ay);

	a *= angularscale/relativeangularwidth2;
	*alpha = a;
	return true;
}

double SquareLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	if (theta.getX() < -angularwidth/2.0)
		return 0;
	if (theta.getX() > angularwidth/2.0)
		return 0;
	if (theta.getY() < -angularwidth/2.0)
		return 0;
	if (theta.getY() > angularwidth/2.0)
		return 0;
	return dens;
}

bool SquareLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	Vector2D<double> t = theta/angularscale;
	axx = (1.0/relativeangularwidth2)*(deriv11(t.getX()-rawd2,t.getY()-rawd2)
			                         -deriv11(t.getX()+rawd2,t.getY()-rawd2)
						 -deriv11(t.getX()-rawd2,t.getY()+rawd2)
						 +deriv11(t.getX()+rawd2,t.getY()+rawd2));
	ayy = (1.0/relativeangularwidth2)*(deriv22(t.getX()-rawd2,t.getY()-rawd2)
			                         -deriv22(t.getX()+rawd2,t.getY()-rawd2)
						 -deriv22(t.getX()-rawd2,t.getY()+rawd2)
						 +deriv22(t.getX()+rawd2,t.getY()+rawd2));
	axy = (1.0/relativeangularwidth2)*(deriv12(t.getX()-rawd2,t.getY()-rawd2)
			                         -deriv12(t.getX()+rawd2,t.getY()-rawd2)
						 -deriv12(t.getX()-rawd2,t.getY()+rawd2)
						 +deriv12(t.getX()+rawd2,t.getY()+rawd2));
	return true;
}

bool SquareLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *potentialval) const
{
	Vector2D<double> t = theta/angularscale;
	double p = potential(t.getX()-rawd2,t.getY()-rawd2)-potential(t.getX()+rawd2,t.getY()-rawd2)
		 - potential(t.getX()-rawd2,t.getY()+rawd2)+potential(t.getX()+rawd2,t.getY()+rawd2);

	*potentialval = (D_ds/D_s)*(1.0/relativeangularwidth2)*angularscale2*p;
	return true;
}

bool SquareLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	*pDeflectionScale = angularscale;
	*pPotentialScale = angularscale*angularscale;
	return true;
}

bool SquareLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 3;
	return true;
}

bool SquareLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	pFloatParams[0] = (float)(ABS(angularwidth)/deflectionScale);
	pFloatParams[1] = (float)(1.0/relativeangularwidth2);
	pFloatParams[2] = (float)(1.0/(relativeangularwidth2*potentialScale)*deflectionScale*deflectionScale);
	return true;
}

std::string SquareLens::getCLProgram(std::string &subRoutineName) const
{
	subRoutineName = "clSquareLensProgram";
	std::string program;

	program += "float clSquareLensPart_phi(float x, float y)\n";
	program += "{\n";
	program += "	float x2 = x*x;\n";
	program += "	float y2 = y*y;\n";
	program += "	return 0.5f*(x2*atan(y/x)+y2*atan(x/y)+x*y*log(x2+y2));\n";
	program += "}\n";
	program += "\n";
	program += "float clSquareLensPart_phix(float x, float y)\n";
	program += "{\n";
	program += "	return x*atan(y/x)+0.5f*y*log(x*x+y*y);\n";
	program += "}\n";
	program += "\n";
	program += "float clSquareLensPart_phiy(float x, float y)\n";
	program += "{\n";
	program += "	return clSquareLensPart_phix(y, x);\n";
	program += "}\n";
	program += "\n";
	program += "float clSquareLensPart_phixx(float x, float y)\n";
	program += "{\n";
	program += "	return atan(y/x);\n";
	program += "}\n";
	program += "\n";
	program += "float clSquareLensPart_phiyy(float x, float y)\n";
	program += "{\n";
	program += "	return atan(x/y);\n";
	program += "}\n";
	program += "\n";
	program += "float clSquareLensPart_phixy(float x, float y)\n";
	program += "{\n";
	program += "	return 0.5f*log(x*x+y*y);\n";
	program += "}\n";
	program += "\n";
	program += "LensQuantities clSquareLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)\n";
	program += "{\n";
	program += "	LensQuantities r;\n";
	program += "	float x = coord.x;\n";
	program += "	float y = coord.y;\n";
	program += "	float a2 = pFloatParams[0]/2.0f;\n";
	program += "	float f = pFloatParams[1];\n";
	program += "	float g = pFloatParams[2];\n";
	program += "\n";
	program += "	r.alphaX = f*(clSquareLensPart_phix(x-a2, y-a2) + clSquareLensPart_phix(x+a2, y+a2) - clSquareLensPart_phix(x+a2, y-a2) - clSquareLensPart_phix(x-a2, y+a2));\n";
	program += "	r.alphaY = f*(clSquareLensPart_phiy(x-a2, y-a2) + clSquareLensPart_phiy(x+a2, y+a2) - clSquareLensPart_phiy(x+a2, y-a2) - clSquareLensPart_phiy(x-a2, y+a2));\n";
	program += "	r.potential = g*(clSquareLensPart_phi(x-a2, y-a2) + clSquareLensPart_phi(x+a2, y+a2) - clSquareLensPart_phi(x+a2, y-a2) - clSquareLensPart_phi(x-a2, y+a2));\n";
	program += "	r.axx = f*(clSquareLensPart_phixx(x-a2, y-a2) + clSquareLensPart_phixx(x+a2, y+a2) - clSquareLensPart_phixx(x+a2, y-a2) - clSquareLensPart_phixx(x-a2, y+a2));\n";
	program += "	r.ayy = f*(clSquareLensPart_phiyy(x-a2, y-a2) + clSquareLensPart_phiyy(x+a2, y+a2) - clSquareLensPart_phiyy(x+a2, y-a2) - clSquareLensPart_phiyy(x-a2, y+a2));\n";
	program += "	r.axy = f*(clSquareLensPart_phixy(x-a2, y-a2) + clSquareLensPart_phixy(x+a2, y+a2) - clSquareLensPart_phixy(x+a2, y-a2) - clSquareLensPart_phixy(x-a2, y+a2));\n";
	program += "\n";
	program += "	return r;\n";
	program += "}\n";

	return program;
}

} // end namespace

