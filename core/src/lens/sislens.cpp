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
#include "sislens.h"
#include "constants.h"
#include "utils.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> SISLensParams::createCopy() const
{
	return std::make_unique<SISLensParams>(dispersion);
}

bool SISLensParams::write(serut::SerializationInterface &si) const
{ 
	if (!si.writeDouble(dispersion)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool SISLensParams::read(serut::SerializationInterface &si)
{ 
	if (!si.readDouble(&dispersion)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

SISLens::SISLens() : SymmetricLens(GravitationalLens::SIS)
{
}

SISLens::~SISLens()
{
}

bool SISLens::processParameters(const GravitationalLensParams *pParams)
{
	const SISLensParams *p = dynamic_cast<const SISLensParams *>(pParams);
	if (!p)
	{
		setErrorString("Parameters are not of type 'SISLensParams'");
		return false;
	}
	
	double sigma = p->getVelocityDispersion();
	
	m_massFactor = CONST_PI*sigma*sigma*getLensDistance()/CONST_G;
	m_densFactor = sigma*sigma/(2.0*CONST_G*getLensDistance());
	m_einstRad = (4.0*CONST_PI/(SPEED_C*SPEED_C)) * sigma*sigma;
	m_velocityDispersion = sigma;
	
	return true;
}

double SISLens::getMassInside(double thetaLength) const
{
	return m_massFactor*thetaLength;
}

double SISLens::getProfileSurfaceMassDensity(double thetaLength) const
{	
	return m_densFactor/thetaLength;
}

bool SISLens::getSurfaceMassDensityDerivative(double thetaLength, double &deriv) const
{
	deriv = -m_densFactor/(thetaLength*thetaLength);
	return true;
}

bool SISLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
	                                double *pPotentialValue) const
{
	*pPotentialValue = (D_ds/D_s)*m_einstRad*theta.getLength();
	return true;
}

bool SISLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	*pDeflectionScale = m_einstRad;
	*pPotentialScale = m_einstRad*m_einstRad;
	return true;
}

bool SISLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 1;
	return true;
}

static const double velScale = 100000.0; // units of 100km/s

bool SISLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	pFloatParams[0] = (float)(m_velocityDispersion/velScale);
	return true;
}

std::string SISLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	std::string program;
	float einstRadFactor = (float)(((4.0*CONST_PI/(SPEED_C*SPEED_C)) * velScale*velScale)/deflectionScale);
	float potentialFactor = (float)(deflectionScale*deflectionScale/potentialScale);

	program += R"XYZ(
LensQuantities clSISLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
	float l = length(coord);
	float l2 = coord.x*coord.x+coord.y*coord.y;
	float scaledVel = pFloatParams[0];
	float scaledEinsteinRad = ((float))XYZ" + float_to_string(einstRadFactor) + R"XYZ()*scaledVel*scaledVel;
	float factor = scaledEinsteinRad/l;
	
	LensQuantities r;
	r.alphaX = factor*coord.x;
	r.alphaY = factor*coord.y;
)XYZ";
	if (potential)
	{
		program += R"XYZ(
	float potentialFactor = (float))XYZ" + float_to_string(potentialFactor) + R"XYZ(;
	r.potential = scaledEinsteinRad*l*potentialFactor;
)XYZ";
	}
	if (derivatives)
	{
		program += R"XYZ(
	r.axx = factor*(1.0f-coord.x*coord.x/l2);
	r.ayy = factor*(1.0f-coord.y*coord.y/l2);
	r.axy = -factor*coord.x*coord.y/l2;
)XYZ";
	}
	program += R"XYZ(
	return r;
}
)XYZ";

	subRoutineName = "clSISLensProgram";
	return program;
}


} // end namespace

