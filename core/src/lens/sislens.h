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
 * \file sislens.h
 */

#ifndef GRALE_SISLENS_H

#define GRALE_SISLENS_H

#include "graleconfig.h"
#include "symmetriclens.h"

namespace grale
{

/** Parameters for a singular isothermal sphere (SIS). */
class GRALE_IMPORTEXPORT SISLensParams : public GravitationalLensParams
{
public:
	SISLensParams()								{ dispersion = 0; }
	
	/** Initializes parameters for a SIS-type lens.
	 *  Initializes parameters for a SIS-type lens.
	 *  \param dispersion Velocity dispersion.
	 */
	SISLensParams(double velocityDispersion) 				{ SISLensParams::dispersion = velocityDispersion; }

	/** Returns the currently set velocity dispersion. */
	double getVelocityDispersion() const					{ return dispersion; }

	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double dispersion;
};

/** A singular isothermal sphere based lens. 
 *  \f[ \Sigma(\theta) = \frac{\sigma_v^2}{2 G D_d \theta} \f]
 *  \f[ M(\theta) = \frac{\pi \sigma_v^2 D_d}{G}\theta \f]
 */
class GRALE_IMPORTEXPORT SISLens : public SymmetricLens
{
public:
	SISLens();
	~SISLens();

	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
	                           double *pPotentialValue) const;

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	double m_massFactor, m_densFactor;
	double m_einstRad;
};

} // end namespace

#endif // GRALE_SISLENS_H

