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
 * \file pointmasslens.h
 */

#ifndef GRALE_POINTMASSLENS_H

#define GRALE_POINTMASSLENS_H

#include "graleconfig.h"
#include "symmetriclens.h"

namespace grale
{

/** Parameters for a point mass lens. */
class GRALE_IMPORTEXPORT PointmassLensParams : public GravitationalLensParams
{
public:
	PointmassLensParams()								{ lensmass = 0; }

	/** Initialize the parameters.
	 *  Initialize the parameters.
	 *  \param mass Mass of the lens.
	 */
	PointmassLensParams(double mass)						{ lensmass = mass; }

	/** Returns the currently set mass. */
	double getLensMass() const							{ return lensmass; }
	
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double lensmass;
};

/** Describes a point mass lens. */
class GRALE_IMPORTEXPORT PointmassLens : public SymmetricLens
{
public:
	PointmassLens();
	~PointmassLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<PointmassLens>(); }
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
	bool getSurfaceMassDensityDerivative(double thetaLength, double &deriv) const override;
private:
	double mass;
};

} // end namespace

#endif // GRALE_POINTMASSLENS_H

