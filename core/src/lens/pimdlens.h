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

#ifndef GRALE_PIMDLENS_H

#define GRALE_PIMDLENS_H

#include "graleconfig.h"
#include "symmetriclens.h"

namespace grale
{

class GRALE_IMPORTEXPORT PIMDLensParams : public GravitationalLensParams
{
public:
	PIMDLensParams();
	PIMDLensParams(double sigma0, double coreRadius, double angularRadius);
	~PIMDLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getCentralDensity() const							{ return m_sigma0; }
	double getCoreRadius() const								{ return m_coreRadius; }
	double getScaleRadius() const								{ return m_scaleRadius; }
private:
	double m_sigma0;
	double m_coreRadius;
	double m_scaleRadius;
};

class GRALE_IMPORTEXPORT PIMDLens : public SymmetricLens
{
public:
	PIMDLens();
	~PIMDLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<PIMDLens>(); }

	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	bool processParameters(const GravitationalLensParams *pLensParams);

	double m_sigma0;
	double m_coreRadius;
	double m_scaleRadius;

	double m_massFactor, m_densFactor;
	double m_a2, m_s2;
};

} // end namespace

#endif // GRALE_PIMDLENS_H
