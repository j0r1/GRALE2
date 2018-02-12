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

#ifndef GRALE_PIEMDLENS_H

#define GRALE_PIEMDLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

class GRALE_IMPORTEXPORT PIEMDLensParams : public GravitationalLensParams
{
public:
	PIEMDLensParams();
	PIEMDLensParams(double sigma0, double coreRadius, double scaleRadius, double epsilon);
	~PIEMDLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	GravitationalLensParams *createCopy() const;

	double getCentralDensity() const										{ return m_sigma0; }
	double getCoreRadius() const											{ return m_coreRadius; }
	double getScaleRadius() const											{ return m_scaleRadius; }
	double getEpsilon() const												{ return m_epsilon; }
private:
	double m_sigma0;
	double m_coreRadius;
	double m_scaleRadius;
	double m_epsilon;
};

class GRALE_IMPORTEXPORT PIEMDLens : public GravitationalLens
{
public:
	PIEMDLens();
	~PIEMDLens();

	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	//bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	Vector2Dd calcI(double omega, Vector2Dd theta) const;
	void getIDerivs(double omega, Vector2D<double> theta, double &axx, double &ayy, double &axy) const;

	double m_sigma0;
	double m_coreRadius;
	double m_scaleRadius;
	double m_epsilon;
	double m_sqrtEpsilon;
	double m_epsFrac;

	double m_angularScale;
	double m_deflectionFactor;
	double m_densFactor;
	double m_a2, m_s2;
};

} // end namespace

#endif // GRALE_PIEMDLENS_H
