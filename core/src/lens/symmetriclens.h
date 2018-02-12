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
 * \file symmetriclens.h
 */

#ifndef GRALE_SYMMETRICLENS_H

#define GRALE_SYMMETRICLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

/** Base class for circularly symmetric lenses. 
 *  \f[ M(\theta) = 2 \pi D_d^2 \int\limits_0^\theta \theta'\Sigma(\theta')\mathrm{d}\theta' \f]
 *  \f[ \vec{\hat{\alpha}}(\vec{\theta}) = \frac{4 G M(\theta)}{c^2 D_d \theta^2}\vec{\theta} \f]
 *  \f[ \hat{\alpha}_{x,x} = \frac{4 G M(\theta)}{c^2 D_d} \frac{\theta_y^2-\theta_x^2}{\theta^4}
 *                         + \frac{8 \pi G D_d \Sigma(\theta)}{c^2}\frac{\theta_x^2}{\theta^2} \f]
 *  \f[ \hat{\alpha}_{y,y} = \frac{4 G M(\theta)}{c^2 D_d} \frac{\theta_x^2-\theta_y^2}{\theta^4}
 *                         + \frac{8 \pi G D_d \Sigma(\theta)}{c^2}\frac{\theta_y^2}{\theta^2} \f]
 *  \f[ \hat{\alpha}_{x,y} = -\frac{4 G M(\theta)}{c^2 D_d} \frac{2\theta_x\theta_y}{\theta^4}
 *                         + \frac{8 \pi G D_d \Sigma(\theta)}{c^2}\frac{\theta_x\theta_y}{\theta^2} \f]
 */
class GRALE_IMPORTEXPORT SymmetricLens : public GravitationalLens
{
protected:
	/** Constructor meant to be used by subclasses. */
	SymmetricLens(GravitationalLens::LensType t);
public:
	~SymmetricLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const			{ return getProfileSurfaceMassDensity(theta.getLength()); }

    static SymmetricLens *cast(GravitationalLens *pLens);
    static const SymmetricLens *cast(const GravitationalLens *pLens);

	/** Implementations need to provide this function.
	 *  Implementations need to provide this function. The function should
	 *  return the total mass within a radius described by \c thetaLength.
	 */
	virtual double getMassInside(double thetaLength) const = 0;
	
	/** Implementations need to provide this function.
	 *  Implementations need to provide this function. The function should
	 *  return the mass density at the radius described by \c thetaLength.
	 */
	virtual double getProfileSurfaceMassDensity(double thetaLength) const = 0;
};

} // end namespace

#endif // GRALE_SYMMETRICLENS_H

