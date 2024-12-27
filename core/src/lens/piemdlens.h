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

#ifndef GRALE_PIEMDLENS_H

#define GRALE_PIEMDLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

class GRALE_IMPORTEXPORT PIEMDLensParams : public GravitationalLensParams
{
public:
	PIEMDLensParams();
	PIEMDLensParams(double sigma0, double coreRadius, double scaleRadius, double epsilon);
	~PIEMDLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getCentralDensity() const										{ return m_sigma0; }
	double getCoreRadius() const											{ return m_coreRadius; }
	double getScaleRadius() const											{ return m_scaleRadius; }
	double getEpsilon() const												{ return m_epsilon; }
private:
	double m_sigma0;
	double m_coreRadius;
	double m_scaleRadius;
	double m_epsilon;
};

class GRALE_IMPORTEXPORT LTPIEMDLensParams : public GravitationalLensParams
{
public:
	LTPIEMDLensParams();
	LTPIEMDLensParams(double velDisp, double coreRadius, double scaleRadius, double ellipticity);
	~LTPIEMDLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getVelocityDispersion() const									{ return m_velDisp; }
	double getCoreRadius() const											{ return m_coreRadius; }
	double getScaleRadius() const											{ return m_scaleRadius; }
	double getEllipticity() const											{ return m_ellipticity; }
private:
	double m_velDisp;
	double m_coreRadius;
	double m_scaleRadius;
	double m_ellipticity;
};

class GRALE_IMPORTEXPORT PIEMDLens : public GravitationalLens
{
public:
	PIEMDLens();
	~PIEMDLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<PIEMDLens>(); }

	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	//bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
	std::unique_ptr<GravitationalLensParams> createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const override;
	std::vector<CLFloatParamInfo> getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const override;
protected:
	PIEMDLens(LensType t);
	bool subInit(double sigma0, double coreRad, double scaleRad, double eps);
	bool processParameters(const GravitationalLensParams *pLensParams);

	Vector2Dd calcI(double omega, Vector2Dd theta) const;
	void getIDerivs(double omega, Vector2D<double> theta, double &axx, double &ayy, double &axy) const;

	double m_coreRadius;
	double m_scaleRadius;
private:
	double m_sigma0;
	double m_epsilon;
	double m_sqrtEpsilon;
	double m_epsFrac;

	double m_angularScale;
	double m_deflectionFactor;
	double m_densFactor;
	double m_a2, m_s2;
};

class GRALE_IMPORTEXPORT LTPIEMDLens : public PIEMDLens
{
public:
	LTPIEMDLens();
	~LTPIEMDLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<LTPIEMDLens>(); }

	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
	std::unique_ptr<GravitationalLensParams> createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const override;
	std::vector<CLFloatParamInfo> getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const override;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double m_velDisp = 0;
	double m_epsHat = 0;
};

} // end namespace

#endif // GRALE_PIEMDLENS_H
