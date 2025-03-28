#pragma once

#include "graleconfig.h"
#include "nfwlens.h"
#include "circularlensprofile.h"

namespace grale
{

class GRALE_IMPORTEXPORT HernquistLensParams : public GravitationalLensParams
{
public:
	HernquistLensParams();
	HernquistLensParams(double sigma_s, double theta_s); // sigma_s is density at theta_s
	~HernquistLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getDensityScale() const										{ return m_densityScale; }
	double getAngularRadiusScale() const								{ return m_angularRadiusScale; }
private:
	double m_densityScale = 0;
	double m_angularRadiusScale = 0;
};

class CircularHernquistLensProfile : public CircularLensProfile
{
public:
	CircularHernquistLensProfile();
	CircularHernquistLensProfile(double sigma_s, double theta_s, double Dd);
	~CircularHernquistLensProfile();

	double getMassInside(double theta) const;
	double getSurfaceMassDensity(double theta) const;
	double getSurfaceMassDensityDerivativeOverTheta(double theta) const;
private:
	static double F(double x) { return NFWLens::F(x); }

	double m_Sigma0 = 0;
	double m_Sigma_s = 0;
	double m_angularRadiusScale = 0;
	double m_massScale = 0;
};

class GRALE_IMPORTEXPORT HernquistLens : public SymmetricLens
{
public:
	HernquistLens();
	~HernquistLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<HernquistLens>(); }

	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	bool processParameters(const GravitationalLensParams *params);

	CircularHernquistLensProfile m_profile;
};

} // end namespace


