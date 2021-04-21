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

#ifndef GRALE_ELLIPTICSERSICLENS_H

#define GRALE_ELLIPTICSERSICLENS_H

#include "graleconfig.h"
#include "ellipticlens.h"

namespace grale
{

class GRALE_IMPORTEXPORT EllipticSersicLensParams : public GravitationalLensParams
{
public:
	EllipticSersicLensParams();
	EllipticSersicLensParams(double centalDensity, double angularScale, double index, double q);
	~EllipticSersicLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getCentralDensity() const							{ return m_centralDensity; }
	double getAngularScale() const								{ return m_angularScale; }
	double getSersicIndex() const								{ return m_sersicIndex; }
	double getEllipticity() const								{ return m_ellipticity; }
private:
	double m_centralDensity;
	double m_angularScale;
	double m_sersicIndex;
	double m_ellipticity;
};

class GRALE_IMPORTEXPORT EllipticSersicLens : public EllipticLens
{
public:
	EllipticSersicLens();
	~EllipticSersicLens();
private:
	bool processParameters(const GravitationalLensParams *pLensParams);

	std::unique_ptr<CircularLensProfile> m_pProfile;
};

} // end namespace

#endif // GRALE_ELLIPTICSERSICLENS_H
