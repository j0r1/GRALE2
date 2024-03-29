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
#include "masssheetlens.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> MassSheetLensParams::createCopy() const
{ 
	return std::make_unique<MassSheetLensParams>(m_density); 
}

bool MassSheetLensParams::write(serut::SerializationInterface &si) const
{ 
	if (!si.writeDouble(m_density)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	return true; 
}

bool MassSheetLensParams::read(serut::SerializationInterface &si)
{ 
	if (!si.readDouble(&m_density)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	return true; 
}

MassSheetLens::MassSheetLens() : GravitationalLens(GravitationalLens::MassSheet)
{
}

MassSheetLens::~MassSheetLens()
{
}

bool MassSheetLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const MassSheetLensParams *pParams = dynamic_cast<const MassSheetLensParams *>(pLensParams);

	if (!pParams)
	{
		setErrorString("Parameters are not of type 'MassSheetLensParams'");
		return false;
	}

	m_density = pParams->getDensity();
	m_factor = (4.0*CONST_PI*CONST_G*m_density*getLensDistance())/(SPEED_C*SPEED_C);
	return true;
}

bool MassSheetLens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
	*pAlpha = theta*m_factor;
	return true;
}

double MassSheetLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	return m_density;
}

bool MassSheetLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	axx = m_factor;
	ayy = m_factor;
	axy = 0;
	return true;
}
	
bool MassSheetLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	*pPotentialValue = 0.5*m_factor*(D_ds/D_s)*theta.getLengthSquared();
	return true;
}

bool MassSheetLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	double angle = ANGLE_ARCSEC;

	*pDeflectionScale = angle;
	*pPotentialScale = angle*angle;
	return true;
}

bool MassSheetLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 2;
	return true;
}

bool MassSheetLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	pFloatParams[0] = (float)(m_factor);
	pFloatParams[1] = (float)(0.5*m_factor*deflectionScale*deflectionScale/potentialScale);
	return true;
}

std::string MassSheetLens::getCLProgram(std::string &subRoutineName, bool derivatives, bool potential) const
{
	std::string program;

	program += "LensQuantities clMassSheetLensProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)\n";
	program += "{\n";
	program += "	float factor = pFloatParams[0];\n";
	program += "	float factor2 = pFloatParams[1];\n";
	program += "\n";
	program += "	LensQuantities r;\n";
	program += "	r.alphaX = factor*coord.x;\n";
	program += "	r.alphaY = factor*coord.y;\n";
	if (potential)
		program += "	r.potential = pFloatParams[1]*(coord.x*coord.x+coord.y*coord.y);\n";
	if (derivatives)
	{
		program += "	r.axx = factor;\n";
		program += "	r.ayy = factor;\n";
		program += "	r.axy = 0;\n";
	}
	program += "	return r;\n";
	program += "}\n";

	subRoutineName = "clMassSheetLensProgram";

	return program;
}

} // end namespace

