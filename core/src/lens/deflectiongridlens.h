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
 * \file deflectiongridlens.h
 */

#ifndef DEFLECTIONGRIDLENS_H

#define DEFLECTIONGRIDLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include <vector>

namespace grale
{

class GridFunction;

class GRALE_IMPORTEXPORT DeflectionGridLensParams : public GravitationalLensParams
{
public:
	DeflectionGridLensParams();
	DeflectionGridLensParams(const std::vector<double> &alphaX, const std::vector<double> &alphaY, 
			         int width, int height, Vector2D<double> bottomLeft, Vector2D<double> topRight);

	const std::vector<double> &getAlphaX() const						{ return m_alphaX; }
	const std::vector<double> &getAlphaY() const						{ return m_alphaY; }
	int getWidth() const									{ return m_width; }
	int getHeight() const									{ return m_height; }
	Vector2D<double> getBottomLeft() const							{ return m_bottomLeft; }
	Vector2D<double> getTopRight() const							{ return m_topRight; }

	GravitationalLensParams *createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	std::vector<double> m_alphaX, m_alphaY;
	int m_width, m_height;
	Vector2D<double> m_topRight, m_bottomLeft;
};

class GRALE_IMPORTEXPORT DeflectionGridLens : public GravitationalLens
{
public:
	DeflectionGridLens();
	~DeflectionGridLens();
	bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const ;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
// 	TODO
//	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);

	GridFunction *m_pAxFunction, *m_pAyFunction;
	std::vector<double> m_alphaX, m_alphaY;
	double m_pixelWidth, m_pixelHeight;
	double m_densFactor;
};

} // end namespace

#endif // DEFLECTIONGRIDLENS_H
