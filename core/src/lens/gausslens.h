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
 * \file gausslens.h
 */

#ifndef GRALE_GAUSSLENS_H

#define GRALE_GAUSSLENS_H

#include "graleconfig.h"
#include "symmetriclens.h"

namespace grale
{

/** Parameters for a lens with a Gaussian 2D mass density distribution. */
class GRALE_IMPORTEXPORT GaussLensParams : public GravitationalLensParams
{
public:
	GaussLensParams()									{ m_mass = 0; m_angularWidth = 0; }
	
	/** Initializes the parameters.
	 *  Initializes the Gaussian lens parameters.
	 *  \param mass Total mass of the lens.
	 *  \param angularwidth Width of the Gaussian distribution.
	 */
	GaussLensParams(double mass, double angularWidth)					{ m_mass = mass; m_angularWidth = angularWidth; }

	/** Returns the currently set mass. */
	double getMass() const 									{ return m_mass; }

	/** Returns the currently set width. */
	double getAngularWidth() const								{ return m_angularWidth; }

	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double m_mass, m_angularWidth;
};

/** Describes a lens with a Gaussian mass distribution. */
class GRALE_IMPORTEXPORT GaussLens : public SymmetricLens
{
public:
	GaussLens();
	~GaussLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<GaussLens>(); }
protected:
	bool processParameters(const GravitationalLensParams *pParams);
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	double mass,angularwidth,angularwidth2;
};

} // end namespace

#endif // GRALE_GAUSSLENS_H

