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

#ifndef GRALE_NFWLENS_H

#define GRALE_NFWLENS_H

#include "graleconfig.h"
#include "symmetriclens.h"

namespace grale
{

class GRALE_IMPORTEXPORT NFWLensParams : public GravitationalLensParams
{
public:
	NFWLensParams();
	NFWLensParams(double rho_s, double theta_s);
	~NFWLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	GravitationalLensParams *createCopy() const;

	double get3DDensityScale() const								{ return m_densityScale3D; }
	double getAngularRadiusScale() const								{ return m_angularRadiusScale; }
private:
	double m_densityScale3D;
	double m_angularRadiusScale;
};

/** Lensing by a spherically symmetric Navarro-Frenk-White profile (NFW).
 *  Lensing by a spherically symmetric Navarro-Frenk-White profile (NFW).
 *
 *  \f[
 *  F(x) = \left\{
 *      \begin{array}{ll} 
 *  	\frac{1}{\sqrt{1-x^2}}\textrm{atanh}{\sqrt{1-x^2}} & x < 1 \\
 *  	1 & x = 1 \\
 *  	\frac{1}{\sqrt{x^2-1}}\textrm{atan}{\sqrt{x^2-1}} & x > 1 \\
 *  	\end{array} \right.
 *  \f]
 *  \f[ \Rightarrow \frac{d F}{d x} = \frac{1-x^2 F(x)}{x(x^2-1)} \f]
 *  \f[ G(x) = \frac{1-F(x)}{x^2-1} \f]
 *  \f[ H(x) = \textrm{ln}\left(\frac{x}{2}\right) + F(x) \Rightarrow \frac{d H}{d x} = G(x) x \f]
 *  \f[ A(x) = \frac{d H}{d x} \f]
 *  \f[ 
 *  B(x) = \left\{
 *      \begin{array}{ll} 
 *  	\textrm{ln}^2\left(\frac{x}{2}\right) - \textrm{atanh}^2\sqrt{1-x^2} & x < 1 \\
 *  	\textrm{ln}^2 2 & x = 1 \\
 *  	\textrm{ln}^2\left(\frac{x}{2}\right) + \textrm{atan}^2\sqrt{x^2-1} & x > 1 \\
 *  	\end{array} \right.
 *  \f]
 *  \f[ \Rightarrow \frac{1}{2} \frac{d B}{d x} = A(x) \f]
 *  \f[ \Sigma(\theta) = 2 r_s \rho_s G\left(\frac{\theta}{\theta_s}\right) \f]
 *  \f[ \textrm{M}(\theta) = 4\pi r_s^3\rho_s H\left(\frac{\theta}{\theta_s}\right) \f]
 *  \f[ \Rightarrow \vec{\hat{\alpha}}\left(\vec{\theta}\right) = \frac{16 \pi G r_s^2 \rho_s}{c^2} A\left(\frac{\theta}{\theta_s}\right) \frac{\vec{\theta}}{\theta} \f]
 *  \f[ \psi(\theta) = \frac{D_{ds}}{D_s} \frac{8 \pi G r_s^2 \rho_s \theta_s}{c^2} B\left(\frac{\theta}{\theta_s}\right) \f]
 *
 *  \b References
 *  \li <A href="http://adsabs.harvard.edu/abs/2001astro.ph..2341K">Keeton, C., A Catalog of Mass Models for Gravitational Lensing, arXiv:astro-ph/0102341, February 2001.</A>
 *  \li <A href="http://adsabs.harvard.edu/abs/2000ApJ...534...34W">Wright, C., Brainerd, T., Gravitational Lensing by NFW Halos, The Astrophysical Journal, Volume 534, Issue 1, pp. 34-40, May 2000.</A>
 */

class GRALE_IMPORTEXPORT NFWLens : public SymmetricLens
{
public:
	NFWLens();
	~NFWLens();

	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	bool processParameters(const GravitationalLensParams *params);
	static double F(double x);
	static double DF(double x);
	static double P(double x);

	double m_deflectionScale;
	double m_angularRadiusScale;
	double m_potentialScale;
	double m_massScale;
	double m_densityScale;
};

inline double NFWLens::F(double x)
{
	if (x < 1.0)
	{
		double tmp = SQRT(1.0-x*x);
		return ATANH(tmp)/tmp;
	}
	else if (x > 1.0)
	{
		double tmp = SQRT(x*x-1.0);
		return ATAN(tmp)/tmp;
	}
	return 1.0;
}

inline double NFWLens::DF(double x)
{
	double x2 = x*x;

	return (1.0-x2*F(x))/(x*(x2-1.0));
}

inline double NFWLens::P(double x)
{
	if (x < 1.0)
	{
		double l = LN(x/2);
		double at = ATANH(SQRT(1.0-x*x));
		return l*l - at*at;
	}
	else if (x > 1.0)
	{
		double l = LN(x/2);
		double at = ATAN(SQRT(x*x-1.0));
		return l*l + at*at;
	}
	return 0.48045301391820142465;
}

} // end namespace

#endif // GRALE_NFWLENS_H

