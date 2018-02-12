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
 * \file sielens.h
 */

#ifndef GRALE_SIELENS_H

#define GRALE_SIELENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

/** Parameters for a singular isothermal ellipse. */
class GRALE_IMPORTEXPORT SIELensParams : public GravitationalLensParams
{
public:
	SIELensParams()								{ V = 0; F = 0; }

	/** Initialize the parameters.
	 *  Initialize the parameters.
	 *  \param v Velocity dispersion.
	 *  \param f Ellipticity factor.
	 * */
	SIELensParams(double velocityDispersion, double ellipticity)		{ V = velocityDispersion; F = ellipticity; }

	/** Returns the currently set velocity dispersion. */
	double getVelocityDispersion() const					{ return V; }

	/** Returns the currently set ellipticity factor. */
	double getEllipticity() const						{ return F; }
	
	GravitationalLensParams *createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double V,F;
};

/** Describes a singular isothermal ellipse. 
 *  \f[ 
 *      \Sigma(\vec{\theta}) = \frac{\sigma_v^2\sqrt{f}}{2 G D_d\sqrt{\theta_x^2+f^2\theta_y^2}}
 *  \f]
 *  \f[
 *  	\vec{\hat{\alpha}}(\vec{\theta}) = \frac{4 \pi \sigma_v^2}{c^2} \frac{\sqrt{f}}{\sqrt{1-f^2}}
 * 					   \left[
 *  	                                       \mathrm{asinh}\left(\frac{\sqrt{1-f^2}}{f}\frac{\theta_x}{|\vec{\theta}|}\right)\vec{e}_x +
 *  	                                       \mathrm{asin}\left(\sqrt{1-f^2}\frac{\theta_y}{|\vec{\theta}|}\right)\vec{e}_y
 *					   \right]
 *  \f]
 *  
 *  \b References
 *  \li <A href="http://adsabs.harvard.edu/abs/1994A&A...284..285K">Kormann, R., Schneider, P., Bartelmann, M., Isothermal elliptical gravitational lens models. Astronomy and Astrophysics, 284:285-299, April 1994.</A>
 */
class GRALE_IMPORTEXPORT SIELens : public GravitationalLens
{
public:
	SIELens();
	~SIELens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double V,F;
};

} // end namespace

#endif // GRALE_SIELENS_H

