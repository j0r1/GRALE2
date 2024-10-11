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

#ifndef GRALE_ALPHAPOTLENS_H

#define GRALE_ALPHAPOTLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

class GRALE_IMPORTEXPORT AlphaPotLensParams : public GravitationalLensParams
{
public:
	AlphaPotLensParams();
	AlphaPotLensParams(double b, double s, double q, double K2, double alpha);
	~AlphaPotLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getB() const { return m_b; }
	double getS() const { return m_s; }
	double getQ() const { return m_q; }
	double getK2() const { return m_K2; }
	double getAlpha() const { return m_alpha; }
private:
	double m_b, m_s, m_q, m_K2, m_alpha;
};

class GRALE_IMPORTEXPORT AlphaPotLens : public GravitationalLens
{
public:
	AlphaPotLens();
	~AlphaPotLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<AlphaPotLens>(); }

	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double m_b, m_s2, m_q2, m_K2, m_alpha;
};

} // end namespace

#endif // GRALE_ALPHAPOTLENS_H
