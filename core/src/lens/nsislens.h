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
 * \file nsislens.h
 */

#ifndef GRALE_NSISLENS_H

#define GRALE_NSISLENS_H

#include "graleconfig.h"
#include "symmetriclens.h"

namespace grale
{
	
/** Parameters for a non-singular isothermal sphere type lens. */
class GRALE_IMPORTEXPORT NSISLensParams : public GravitationalLensParams
{
public:
	NSISLensParams()							{ dispersion = 0; corerad = 0; }

	/** Initialize the parameters.
	 *  Initialize the parameters.
	 *  \param dispersion Velocity dispersion.
	 *  \param angularcorerad Angular core radius.
	 */
	NSISLensParams(double velocityDispersion, double angularCoreRadius)	{ NSISLensParams::dispersion = velocityDispersion; corerad = angularCoreRadius; }

	/** Returns the currently set velocity dispersion. */
	double getVelocityDispersion() const					{ return dispersion; }

	/** Returns the currently set angular core radius. */
	double getAngularCoreRadius() const					{ return corerad; }

	GravitationalLensParams *createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double dispersion,corerad;
};

/** Describes a non-singular isothermal sphere type lens. 
 *  \f[ \Sigma(\theta) = \frac{\sigma_v^2}{2 G D_d \sqrt{\theta^2 + \theta_c^2}} \f]
 *  \f[ M(\theta) = \frac{\pi \sigma_v^2 D_d}{G}\left(\sqrt{\theta^2+\theta_c^2}-\theta_c\right)\f]
 */
class GRALE_IMPORTEXPORT NSISLens : public SymmetricLens
{
public:
	NSISLens();
	~NSISLens();
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	double m_massFactor, m_densFactor, m_core, m_coreSquared;
};

} // end namespace

#endif // GRALE_NSISLENS_H

