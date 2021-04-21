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

#ifndef GRALE_ELLIPTICLENS_H

#define GRALE_ELLIPTICLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

class Real2DFunction;
class CircularLensProfile;

/** Base class for the lens effect by an elliptical lens, based on the information of the circular equivalent.
 *  Base class for the lens effect by an elliptical lens, based on the information of the circular equivalent.
 *   
 *  Suppose we have a circularly symmetric lens, in which the integrated (projected) mass is described by:
 *  \f[ 
 *  	\textrm{M}_c(\theta) = 2\pi D_d^2 \int_0^\theta \Sigma_c(\theta') \theta' d\theta' 
 *  \f]
 *
 *  The lens effect by the mass distribution
 *  \f[ \Sigma\left(\vec{\theta}\right) = \Sigma_c\left(\sqrt{\theta_x^2+\frac{\theta_y^2}{q^2}}\right) \f]
 *  can then be described as follows.
 *
 *  First, define the following quantities:
 *  \f[ \xi(u,x,y) = \sqrt{u\left(x^2+\frac{y^2}{1-(1-q^2)u}\right)}  \f]
 *  \f[ \textrm{H}(\theta) = \frac{1}{\theta}\frac{d\Sigma_c}{d\theta}(\theta) \f]
 *  \f[	\textrm{I}(x,y) = \frac{4 G}{c^2 D_d}\int_0^1 \frac{\textrm{M}_c\left(\xi(u,x,y)\right)}{u} \frac{1}{\sqrt{1-(1-q^2)u}} du \f]
 *  \f[ \textrm{J}_n(x,y) = \frac{4 \pi G D_d}{c^2} \int_0^1 \frac{\Sigma_c\left(\xi(u,x,y)\right)}{\left(1-(1-q^2)u\right)^{\left(n+\frac{1}{2}\right)}} du \f]
 *  \f[ \textrm{K}_n(x,y) = \frac{4 \pi G D_d}{c^2} \int_0^1 \textrm{H}\left(\xi(u,x,y)\right) \frac{u}{\left(1-(1-q^2)u\right)^{\left(n+\frac{1}{2}\right)}} du \f]
 *
 *  One then has the following expressions for the lensing potential, deflection angle and its derivatives:
 *  \f[ \psi\left(\vec{\theta}\right) = \frac{D_{ds}}{D_s} \frac{q}{2} \textrm{I}(\theta_x, \theta_y) \f]
 *  \f[ \hat{\alpha}_x\left(\vec{\theta}\right) = q \theta_x \textrm{J}_0(\theta_x, \theta_y) \f]
 *  \f[ \hat{\alpha}_y\left(\vec{\theta}\right) = q \theta_y \textrm{J}_1(\theta_x, \theta_y) \f]
 *  \f[ \frac{\partial \hat{\alpha}_x}{\partial\theta_x}\left(\vec{\theta}\right) = q \textrm{J}_0(\theta_x, \theta_y) + q \theta_x^2 \textrm{K}_0(\theta_x, \theta_y) \f]
 *  \f[ \frac{\partial \hat{\alpha}_y}{\partial\theta_y}\left(\vec{\theta}\right) = q \textrm{J}_1(\theta_x, \theta_y) + q \theta_y^2 \textrm{K}_2(\theta_x, \theta_y) \f]
 *  \f[ \frac{\partial \hat{\alpha}_x}{\partial\theta_y}\left(\vec{\theta}\right) = q \theta_x\theta_y \textrm{K}_1(\theta_x, \theta_y) \f]
 *
 *  \b References
 *  \li <A href="http://adsabs.harvard.edu/abs/2001astro.ph..2341K">Keeton, C., A Catalog of Mass Models for Gravitational Lensing, arXiv:astro-ph/0102341, February 2001.</A>
 */
class GRALE_IMPORTEXPORT EllipticLens : public GravitationalLens
{
protected:
	EllipticLens(GravitationalLens::LensType t);
public:
	~EllipticLens();

	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	void subInit(double q, CircularLensProfile *pProfile,
		           double absError,
			   double relError,
			   int limit);
private:
	std::unique_ptr<Real2DFunction> m_pI, m_pJ0, m_pJ1, m_pK0, m_pK1, m_pK2;
	CircularLensProfile *m_pProfile;
	double m_q;
};

} // end namespace

#endif // GRALE_ELLIPTICLENS_H

