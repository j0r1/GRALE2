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

#ifndef GRALE_PROFILELENS_H

#define GRALE_PROFILELENS_H

#include "graleconfig.h"
#include "symmetriclens.h"
#include <vector>

namespace grale
{

class GRALE_IMPORTEXPORT ProfileLensParams : public GravitationalLensParams
{
public:
	ProfileLensParams()									{ }
	ProfileLensParams(double endRadius, const std::vector<double> &profile)			{ m_profile = profile; m_endRadius = endRadius; }
	const std::vector<double> &getProfile() const						{ return m_profile; }
	double getEndRadius() const								{ return m_endRadius; }
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double m_endRadius;
	std::vector<double> m_profile;
};

class GRALE_IMPORTEXPORT ProfileLens : public SymmetricLens
{
public:
	ProfileLens();
	~ProfileLens();
	
	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<ProfileLens>(); }

	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
		                   double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	std::vector<double> m_profile;
	std::vector<double> m_mass;
	std::vector<double> m_potential;
	double m_endRadius;
	double m_stepSize;
};

} // end namespace

#endif // GRALE_PROFILELENS_H
