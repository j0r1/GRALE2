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
 * \file harmoniclens.h
 */

#ifndef GRALE_HARMONICLENS_H

#define GRALE_HARMONICLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

class GRALE_IMPORTEXPORT HarmonicLensParams : public GravitationalLensParams
{
public:
	HarmonicLensParams()																{ m_sigma0 = 0; m_k = 1; m_l = 1; m_phiX = 0; m_phiY = 0; }
	HarmonicLensParams(double sigma0, double k, double l, double phiX = 0, double phiY = 0)	{ m_sigma0 = sigma0; m_k = k; m_l = l; m_phiX = phiX; m_phiY = phiY; }

	double getDensityScale() const														{ return m_sigma0; }
	double getK() const																	{ return m_k; }
	double getL() const																	{ return m_l; }
	double getPhiX() const																{ return m_phiX; }
	double getPhiY() const																{ return m_phiY; }
	
	GravitationalLensParams *createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double m_sigma0, m_phiX, m_phiY, m_factor;
	double m_k, m_l;
};

class GRALE_IMPORTEXPORT HarmonicLens : public GravitationalLens
{
public:
	HarmonicLens();
	~HarmonicLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double m_sigma0, m_phiX, m_phiY, m_factor;
	double m_k, m_l;
};

} // end namespace

#endif // GRALE_HARMONICLENS_H

