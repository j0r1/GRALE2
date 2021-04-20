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
 * \file nsielens.h
 */

#ifndef GRALE_NSIELENS_H

#define GRALE_NSIELENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{
	
/** Parameters for a non-singular isothermal ellipse. */
class GRALE_IMPORTEXPORT NSIELensParams : public GravitationalLensParams
{
public:
	NSIELensParams()							{ m_velocityDispersion = 0; m_ellipticity = 0; m_angularCoreRadius = 0; }
	/** Initialize the parameters.
	 *  Initialize the parameters.
	 *  \param v Velocity dispersion.
	 *  \param f Ellipticity factor.
	 *  \param angularcoreradius Angular core radius.
	 */
	NSIELensParams(double velocityDispersion, double ellipticity, double angularCoreRadius)
										{ m_velocityDispersion = velocityDispersion; m_ellipticity = ellipticity; m_angularCoreRadius = angularCoreRadius; }

	/** Returns the currently set velocity dispersion. */
	double getVelocityDispersion() const					{ return m_velocityDispersion; }

	/** Returns the currently set ellipticity factor. */
	double getEllipticity() const						{ return m_ellipticity; }

	/** Returns the currently set angular core radius. */
	double getAngularCoreRadius() const					{ return m_angularCoreRadius; }
	
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double m_velocityDispersion;
	double m_ellipticity;
	double m_angularCoreRadius;
};

/** Describes a non-singular isothermal ellipse. 
 *  \f[ \Sigma(\vec{\theta}) = \frac{\sigma_v^2 \sqrt{f}}{2 G D_d \sqrt{\theta_x^2+f^2\theta_y^2+\theta_c^2}} \f]
 *  \f[ \hat{\alpha}_x(\vec{\theta}) = \frac{4\pi\sigma_v^2}{c^2}\frac{\sqrt{f}}{\sqrt{1-f^2}}
 *                       \mathrm{atanh}\left(\frac{\theta_x\sqrt{1-f^2}}{\sqrt{\theta_x^2+f^2\theta_y^2+\theta_c^2}+f\theta_c}\right) \f]
 *  \f[ \hat{\alpha}_y(\vec{\theta}) = \frac{4\pi\sigma_v^2}{c^2}\frac{\sqrt{f}}{\sqrt{1-f^2}}
 *                       \mathrm{atan}\left(\frac{f\theta_y\sqrt{1-f^2}}{f\sqrt{\theta_x^2+f^2\theta_y^2+\theta_c^2}+\theta_c}\right) \f]
 *  
 *  \b References
 *     \li <A href="http://adsabs.harvard.edu/abs/2004astro.ph..7232K">Kochanek, C.S., Schneider, P., Wambsganss, J., 2004, Part 2 of Gravitational Lensing: Strong, Weak & Micro, Proceedings of the 33rd Saas-Fee Advanced Course, G. Meylan, P. Jetzer & P. North, eds. (Springer-Verlag: Berlin) </A>
 */
class GRALE_IMPORTEXPORT NSIELens : public GravitationalLens
{
public:
	NSIELens();
	~NSIELens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double V,F,corerad;
};

} // end namespace

#endif // GRALE_NSIELENS_H

