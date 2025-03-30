#pragma once

#include "graleconfig.h"
#include "ellipticlens.h"

namespace grale
{

class GRALE_IMPORTEXPORT EllipticHernquistLensParams : public GravitationalLensParams
{
public:
	EllipticHernquistLensParams();
	EllipticHernquistLensParams(double sigma_s, double theta_s, double q);
	~EllipticHernquistLensParams();

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::unique_ptr<GravitationalLensParams> createCopy() const;

	double getDensityScale() const										{ return m_densityScale; }
	double getAngularRadiusScale() const								{ return m_angularRadiusScale; }
	double getEllipticity() const										{ return m_ellipticity; }
private:
	double m_densityScale = 0;
	double m_angularRadiusScale = 0;
	double m_ellipticity = 0;
};

class GRALE_IMPORTEXPORT EllipticHernquistLens : public EllipticLens
{
public:
	EllipticHernquistLens();
	~EllipticHernquistLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<EllipticHernquistLens>(); }
private:
	bool processParameters(const GravitationalLensParams *pLensParams);

	std::unique_ptr<CircularLensProfile> m_pProfile;
};

} // end namespace

