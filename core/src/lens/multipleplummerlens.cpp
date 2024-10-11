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
#include "multipleplummerlens.h"
#include "constants.h"
#include "utils.h"
#include <iostream>

// TODO: FOR DEBUGGING
#include <stdio.h>

namespace grale
{

std::unique_ptr<GravitationalLensParams> MultiplePlummerLensParams::createCopy() const
{
	return std::make_unique<MultiplePlummerLensParams>(m_lensInfo);
}

bool MultiplePlummerLensParams::write(serut::SerializationInterface &si) const
{
	int32_t numlenses = (int32_t)m_lensInfo.size();
	
	if (!si.writeInt32(numlenses))
	{
		setErrorString(std::string("Error writing number of lenses: ") + si.getErrorString());
		return false;
	}

	for (auto it = m_lensInfo.begin() ; it != m_lensInfo.end() ; it++)
	{
		double array[4];
		
		array[0] = (*it).getMass();
		array[1] = (*it).getAngularWidth();
		array[2] = (*it).getAngularPosition().getX();
		array[3] = (*it).getAngularPosition().getY();
		
		if (!si.writeDoubles(array,4))
		{
			setErrorString(std::string("Error writing plummer lens parameters: ") + si.getErrorString());
			return false;
		}
	}
	return true;
}

bool MultiplePlummerLensParams::read(serut::SerializationInterface &si)
{
	int32_t numlenses,i;

	if (!si.readInt32(&numlenses))
	{
		setErrorString(std::string("Error reading number of lenses: ") + si.getErrorString());
		return false;
	}
	
	m_lensInfo.clear();
	for (i = 0 ; i < numlenses ; i++)
	{
		double array[4];

		if (!si.readDoubles(array,4))
		{
			setErrorString(std::string("Error reading plummer lens parameters: ") + si.getErrorString());
			return false;
		}
		
		Vector2D<double> v(array[2],array[3]);
		PlummerLensInfo p(array[0],array[1],v);
		m_lensInfo.push_back(p);
	}
	
	return true;
}

MultiplePlummerLens::MultiplePlummerLens() : GravitationalLens(GravitationalLens::MultiplePlummers)
{
}

MultiplePlummerLens::~MultiplePlummerLens()
{
}

bool MultiplePlummerLens::processParameters(const GravitationalLensParams *params)
{
	const MultiplePlummerLensParams *p = dynamic_cast<const MultiplePlummerLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'MultiplePlummerLensParams'");
		return false;
	}
	
	int i;

	numlenses = p->getLensInfo().size();
	lensinfo.resize(numlenses);

	totalmass = 0;
	auto it = p->getLensInfo().begin();
	for (i = 0 ; i < numlenses ; i++,it++)
		totalmass += (*it).getMass();

	totalmass = ABS(totalmass);
	scalefactor = SQRT((4.0*CONST_G*totalmass)/(SPEED_C*SPEED_C*getLensDistance()));
		
	for (i = 0, it = p->getLensInfo().begin() ; i < numlenses ; i++,it++)
		lensinfo[i] = PlummerLensInfo((*it).getMass()/totalmass,(*it).getAngularWidth()/scalefactor,(*it).getAngularPosition());
	
	return true;
}

bool MultiplePlummerLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *alpha) const
{	
	Vector2D<double> sum(0,0);
	
	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[i].getAngularWidth();
		double denom = scaleddiff.getLengthSquared()+w*w;

		sum += diff*lensinfo[i].getMass()/denom;
	}

	*alpha = sum;
	return true;
}

double MultiplePlummerLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double sum = 0;
	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[i].getAngularWidth();
		Vector2D<double> s = scaleddiff/w;
		double denom = w*(s.getLengthSquared()+1.0);

		sum += lensinfo[i].getMass()/(denom*denom);
	}
	return (sum*SPEED_C*SPEED_C)/(4.0*CONST_PI*CONST_G*getLensDistance());
}

bool MultiplePlummerLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	axx = getDerivAlphaTheta(1,1,theta);
	ayy = getDerivAlphaTheta(2,2,theta);
	axy = getDerivAlphaTheta(1,2,theta);
	return true;
}

bool MultiplePlummerLens::getAlphaVectorSecondDerivatives(Vector2D<double> theta, double &axxx, double &ayyy, double &axxy, double &ayyx) const
{
	double axxxSum = 0, ayyySum = 0, axxySum = 0, ayyxSum = 0;

	for (int k = 0 ; k < numlenses ; k++)
	{
		Vector2Dd diff = theta-lensinfo[k].getAngularPosition();
		Vector2Dd scaledDiff = diff/scalefactor;
		double w = lensinfo[k].getAngularWidth();
		double denom1 = scaledDiff.getLengthSquared()+w*w;
		double denom2 = denom1*denom1;
		double denom3 = denom2*denom1;
		double theta_x = scaledDiff.getX();
		double theta_y = scaledDiff.getY();
		double theta_x2 = theta_x*theta_x;
		double theta_y2 = theta_y*theta_y;
		double theta_x3 = theta_x2*theta_x;
		double theta_y3 = theta_y2*theta_y;
		double frac = lensinfo[k].getMass();

		axxxSum += frac*(8.0*theta_x3/denom3 - 6.0*theta_x/denom2);
		ayyySum += frac*(8.0*theta_y3/denom3 - 6.0*theta_y/denom2);
		axxySum += frac*(8.0*theta_x2*theta_y/denom3 - 2.0*theta_y/denom2);
		ayyxSum += frac*(8.0*theta_y2*theta_x/denom3 - 2.0*theta_x/denom2);
	}

	axxx = axxxSum/scalefactor;
	ayyy = ayyySum/scalefactor;
	axxy = axxySum/scalefactor;
	ayyx = ayyxSum/scalefactor;
	return true;
}

double MultiplePlummerLens::getDerivAlphaTheta(int i,int j,Vector2D<double> theta) const
{
	double deriv = 0;
	double deltaij = 0;

	if (i == j)
		deltaij = 1;
	
	for (int k = 0 ; k < numlenses ; k++)
	{
		Vector2D<double> diff = theta-lensinfo[k].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[k].getAngularWidth();
		double denom = scaleddiff.getLengthSquared()+w*w;
		double term1 = deltaij/denom;
		double term2 = (2.0*scaleddiff.getComponent(i)*scaleddiff.getComponent(j))/(denom*denom);
		
		deriv += lensinfo[k].getMass()*(term1-term2);
	}
	return deriv;
}

bool MultiplePlummerLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *potentialval) const
{
	double sum = 0;
	
	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[i].getAngularWidth();
		double a = scaleddiff.getLengthSquared()+w*w;
	
		sum += (LN(a)+2.0*LN(scalefactor))*lensinfo[i].getMass();
	}

	*potentialval = (D_ds/D_s)*((2.0*CONST_G*totalmass)/((SPEED_C*SPEED_C)*getLensDistance()))*sum;
	return true;
}

bool MultiplePlummerLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	*pDeflectionScale = SQRT((4.0*CONST_G*totalmass)/(SPEED_C*SPEED_C*getLensDistance()));
	*pPotentialScale = (4.0*CONST_G*totalmass)/(SPEED_C*SPEED_C*getLensDistance());
	return true;
}

bool MultiplePlummerLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 4*numlenses;
	return true;
}

bool MultiplePlummerLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	for (int i = 0 ; i < numlenses ; i++)
	{
		double scaledWidth = lensinfo[i].getAngularWidth()*scalefactor/deflectionScale;

		// mass factor
		pFloatParams[i*4+0] = (float)((4.0*CONST_G*lensinfo[i].getMass()*totalmass)/(SPEED_C*SPEED_C*getLensDistance()*deflectionScale*deflectionScale));
		// width factor
		pFloatParams[i*4+1] = (float)(scaledWidth);
		// scaled position
		pFloatParams[i*4+2] = (float)(lensinfo[i].getAngularPosition().getX()/deflectionScale);
		pFloatParams[i*4+3] = (float)(lensinfo[i].getAngularPosition().getY()/deflectionScale);
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> MultiplePlummerLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	auto params = std::make_unique<MultiplePlummerLensParams>();
	for (int i = 0 ; i < numlenses ; i++)
	{
		float massparam = pFloatParams[i*4+0];
		float widthparam = pFloatParams[i*4+1];
		float x = pFloatParams[i*4+2];
		float y = pFloatParams[i*4+3];

		double mass = ((double)massparam)*(SPEED_C*SPEED_C*getLensDistance()*deflectionScale*deflectionScale)/(4.0*CONST_G);
		double width = widthparam*deflectionScale;
		Vector2Dd pos(deflectionScale*(double)x, deflectionScale*(double)y);
		params->addLensInfo({mass, width, pos});
	}

	return params;
}

std::string MultiplePlummerLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	std::string prog;
	float potentialPrefactor = (float)(deflectionScale*deflectionScale/(2.0*potentialScale));

	prog += R"XYZ(
LensQuantities clMultiplePlummerLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
	LensQuantities r = { 0 } ;

	for (int i = 0 ; i < )XYZ" + std::to_string(numlenses) +  R"XYZ( ; i++)
	{
		int i4 = i*4;
		float scaledMass = pFloatParams[i4+0];
		float scaledWidth = pFloatParams[i4+1];
		float posx = pFloatParams[i4+2];
		float posy = pFloatParams[i4+3];
		float dx = coord.x-posx;
		float dy = coord.y-posy;
		float dx2 = dx*dx;
		float dy2 = dy*dy;
		float w2 = scaledWidth*scaledWidth;
		float denom = dx2+dy2+w2;
		float denom2 = denom*denom;
		float factor = scaledMass/denom;
		float factor2 = scaledMass/denom2;
	
		r.alphaX += factor*dx;
		r.alphaY += factor*dy;
)XYZ";
	if (potential)
		prog += R"XYZ(
		float potentialPrefactor = (float))XYZ" + float_to_string(potentialPrefactor) + R"XYZ(;
		r.potential += potentialPrefactor*scaledMass*log(denom);
)XYZ";
	if (derivatives)
	{
		prog += R"XYZ(
		r.axx += factor2*(-dx2+dy2+w2);
		r.ayy += factor2*(+dx2-dy2+w2);
		r.axy += factor2*(-2.0*dx*dy);
)XYZ";
	}
	prog += R"XYZ(
	}
	return r;
}
)XYZ";

	subRoutineName = "clMultiplePlummerLensProgram";

	return prog;
}

} // end namespace

