#include "graleconfig.h"
#include "elliptichernquistlens.h"
#include "constants.h"
#include "hernquistlens.h"

namespace grale
{

EllipticHernquistLensParams::EllipticHernquistLensParams()
{
}

EllipticHernquistLensParams::EllipticHernquistLensParams(double sigma_s, double theta_s, double q)
	: m_densityScale(sigma_s),
	  m_angularRadiusScale(theta_s),
	  m_ellipticity(q)
{
}

EllipticHernquistLensParams::~EllipticHernquistLensParams()
{
}

bool EllipticHernquistLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_densityScale) && si.writeDouble(m_angularRadiusScale) &&
	      si.writeDouble(m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool EllipticHernquistLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_densityScale) && si.readDouble(&m_angularRadiusScale) &&
	      si.readDouble(&m_ellipticity)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> EllipticHernquistLensParams::createCopy() const
{
	return std::make_unique<EllipticHernquistLensParams>(m_densityScale, m_angularRadiusScale, m_ellipticity);
}

EllipticHernquistLens::EllipticHernquistLens() : EllipticLens(EllipticHernquist)
{
}

EllipticHernquistLens::~EllipticHernquistLens()
{
}

bool EllipticHernquistLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const EllipticHernquistLensParams *pParams = dynamic_cast<const EllipticHernquistLensParams *>(pLensParams);

	if (pParams == 0)
	{
		setErrorString("Lens parameters are not of type 'EllipticHernquistLensParams'");
		return false;
	}

	double sigma_s = pParams->getDensityScale();
	double theta_s = pParams->getAngularRadiusScale();
	double q = pParams->getEllipticity();
	double D_d = getLensDistance();

	m_pProfile = std::make_unique<CircularHernquistLensProfile>(sigma_s, theta_s, D_d);

	subInit(q, m_pProfile.get(), 0, 1e-5, 1024);

	return true;
}

bool EllipticHernquistLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	setErrorString("Deflection derivatives currently not supported (you can create a deflection grid lens if needed)");
	return false;
}

} // end namespace
