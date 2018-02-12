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

#ifndef GRALE_ELLIPTICNFWLENS_H

#define GRALE_ELLIPTICNFWLENS_H

#include "graleconfig.h"
#include "ellipticlens.h"

namespace grale
{

class GRALE_IMPORTEXPORT EllipticNFWLensParams : public GravitationalLensParams
{
public:
	EllipticNFWLensParams();
	EllipticNFWLensParams(double rho_s, double theta_s, double q);
	~EllipticNFWLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	GravitationalLensParams *createCopy() const;

	double get3DDensityScale() const								{ return m_densityScale3D; }
	double getAngularRadiusScale() const								{ return m_angularRadiusScale; }
	double getEllipticity() const									{ return m_ellipticity; }
private:
	double m_densityScale3D;
	double m_angularRadiusScale;
	double m_ellipticity;
};

class GRALE_IMPORTEXPORT EllipticNFWLens : public EllipticLens
{
public:
	EllipticNFWLens();
	~EllipticNFWLens();
private:
	bool processParameters(const GravitationalLensParams *pLensParams);

	CircularLensProfile *m_pProfile;
};

} // end namespace

#endif // GRALE_ELLIPTICNFWLENS_H
