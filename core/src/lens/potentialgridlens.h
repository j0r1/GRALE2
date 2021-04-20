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
 * \file potentialgridlens.h
 */

#ifndef GRALE_POTENTIALGRIDLENS_H

#define GRALE_POTENTIALGRIDLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include <gsl/gsl_interp2d.h>

namespace grale
{

class GRALE_IMPORTEXPORT PotentialGridLensParams : public GravitationalLensParams
{
public:
	PotentialGridLensParams();
	PotentialGridLensParams(Vector2Dd bottomLeft, Vector2Dd topRight, const std::vector<double> &values,
			int numX, int numY);
	~PotentialGridLensParams();

	Vector2Dd getBottomLeft() const								{ return m_bottomLeft; }
	Vector2Dd getTopRight() const								{ return m_topRight; }
	const std::vector<double> &getValues() const				{ return m_values; }
	int getNumX() const											{ return m_numX; }
	int getNumY() const											{ return m_numY; }

	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	Vector2Dd m_bottomLeft, m_topRight;
	std::vector<double> m_values;
	int m_numX, m_numY;
};

class GRALE_IMPORTEXPORT PotentialGridLens : public GravitationalLens
{
public:
	PotentialGridLens();
	~PotentialGridLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	gsl_interp2d *m_pInterp;
	gsl_interp_accel *m_pXAccel, *m_pYAccel;

	std::vector<double> m_x, m_y, m_z;
};

} // end namespace

#endif // GRALE_POTENTIALGRIDLENS_H

