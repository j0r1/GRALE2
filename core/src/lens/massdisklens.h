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
 * \file massdisklens.h
 */

#ifndef GRALE_MASSDISKLENS_H

#define GRALE_MASSDISKLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include "constants.h"

namespace grale
{

class GRALE_IMPORTEXPORT MassDiskLensParams : public GravitationalLensParams
{
public:
	MassDiskLensParams()								{ m_density = 0; m_angularRadius = 0; }
	MassDiskLensParams(double density, double angularRadius)			{ m_density = density; m_angularRadius = angularRadius; }
	MassDiskLensParams(double Dd, double Ds, double Dds, double angularRadius)	{ m_density = (SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*Dd))*(Ds/Dds); m_angularRadius = angularRadius; }
	double getDensity() const							{ return m_density; }
	double getAngularRadius() const							{ return m_angularRadius; }
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double m_density, m_angularRadius; 
};

class GRALE_IMPORTEXPORT MassDiskLens : public GravitationalLens
{
public:
	MassDiskLens();
	~MassDiskLens();

	bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;

	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
	                           double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double m_factor;
	double m_scaleFactor;
	double m_density;
	double m_angularRadius; 
	double m_angularRadiusSquared; 
	double m_scaledMass;
};

} // end namespace

#endif // GRALE_MASSDISKLENS_H
