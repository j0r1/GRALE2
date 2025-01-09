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

class GRALE_IMPORTEXPORT LTPIMDLensParams : public GravitationalLensParams
{
public:
	LTPIMDLensParams();
	LTPIMDLensParams(double velDisp, double coreRadius, double angularRadius);
	~LTPIMDLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getVelocityDispersion() const						{ return m_velDisp; }
	double getCoreRadius() const								{ return m_coreRadius; }
	double getScaleRadius() const								{ return m_scaleRadius; }
private:
	double m_velDisp;
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

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const override;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const override;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const override;
	std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
	std::unique_ptr<GravitationalLensParams> createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const override;
	std::vector<CLFloatParamInfo> getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const override;
protected:
	PIMDLens(LensType t);
	bool subInit(double sigma0, double coreRad, double scaleRad);
	bool processParameters(const GravitationalLensParams *pLensParams);
	std::string getPIMDMainProgram(const std::string &subRoutineName, const std::string &pimdParams, double deflectionScale, double potentialScale, bool derivatives, bool potential) const;

	double m_sigma0;
	double m_coreRadius;
	double m_scaleRadius;

	double m_massFactor, m_densFactor;
	double m_a2, m_s2;
	double m_angularScale;
};

class GRALE_IMPORTEXPORT LTPIMDLens : public PIMDLens
{
public:
	LTPIMDLens();
	~LTPIMDLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<LTPIMDLens>(); }

	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const override;
	std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
	std::unique_ptr<GravitationalLensParams> createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const override;
	std::vector<CLFloatParamInfo> getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const override;
private:
	bool processParameters(const GravitationalLensParams *pLensParams);

	double m_velDisp = 0;
};

} // end namespace

#endif // GRALE_PIMDLENS_H
