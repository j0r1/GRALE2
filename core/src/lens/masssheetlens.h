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

/**
 * \file masssheetlens.h
 */

#ifndef GRALE_MASSSHEETLENS_H

#define GRALE_MASSSHEETLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include "constants.h"

namespace grale
{

class GRALE_IMPORTEXPORT MassSheetLensParams : public GravitationalLensParams
{
public:
	MassSheetLensParams()								{ m_density = 0; }
	MassSheetLensParams(double density)						{ m_density = density; }
	MassSheetLensParams(double Dd, double Ds, double Dds)				{ m_density = (SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*Dd))*(Ds/Dds); }
	double getDensity() const							{ return m_density; }
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double m_density; 
};

class GRALE_IMPORTEXPORT MassSheetLens : public GravitationalLens
{
public:
	MassSheetLens();
	~MassSheetLens();

	bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;

	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
	                           double *pPotentialValue) const;

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double m_factor;
	double m_density;
};

} // end namespace

#endif // GRALE_MASSSHEETLENS_H
