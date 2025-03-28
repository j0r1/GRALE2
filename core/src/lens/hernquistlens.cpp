#include "graleconfig.h"
#include "hernquistlens.h"
#include "constants.h"

#include <iostream>

namespace grale
{

HernquistLensParams::HernquistLensParams()
{
}

HernquistLensParams::HernquistLensParams(double sigma_s, double theta_s)
{
	m_densityScale = sigma_s;
	m_angularRadiusScale = theta_s;
}

HernquistLensParams::~HernquistLensParams()
{
}

bool HernquistLensParams::write(serut::SerializationInterface &si) const
{
	if (!(si.writeDouble(m_densityScale) && si.writeDouble(m_angularRadiusScale)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool HernquistLensParams::read(serut::SerializationInterface &si)
{
	if (!(si.readDouble(&m_densityScale) && si.readDouble(&m_angularRadiusScale)))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

std::unique_ptr<GravitationalLensParams> HernquistLensParams::createCopy() const
{
	return std::make_unique<HernquistLensParams>(m_densityScale, m_angularRadiusScale);
}

HernquistLens::HernquistLens() : SymmetricLens(GravitationalLens::Hernquist)
{
}

HernquistLens::~HernquistLens()
{
}

bool HernquistLens::processParameters(const GravitationalLensParams *p)
{
	const HernquistLensParams *pParams = dynamic_cast<const HernquistLensParams *>(p);
	if (!pParams)
	{
		setErrorString("Parameters are not of type 'HernquistLensParams'");
		return false;
	}

	double sigma_s = pParams->getDensityScale();
	double theta_s = pParams->getAngularRadiusScale();
	double r_s = theta_s * getLensDistance();
	double Dd = getLensDistance();

	m_angularRadiusScale = theta_s;
	m_Sigma_s = sigma_s;
	m_Sigma0 = 15.0*sigma_s/4.0;
	m_massScale = 2.0*CONST_PI*Dd*Dd*theta_s*theta_s*m_Sigma0;

	return true;
}

double HernquistLens::getMassInside(double thetaLength) const
{
	double x = thetaLength/m_angularRadiusScale;
	double x2 = x*x;
	double gx = -2.0/3.0;

	if (x < 1.0)
	{
		double oneMinX2 = 1.0-x2;
		gx = 1.0/(x2-1.0) + x2*std::atanh(std::sqrt(oneMinX2))/std::pow(oneMinX2, 1.5);
	}
	else if (x > 1.0)
	{
		double x2MinOne = x2-1.0;
		gx = 1.0/x2MinOne - x2*std::atan(std::sqrt(x2MinOne))/std::pow(x2MinOne, 1.5);
	}

	return m_massScale * (gx + 1.0);
}

double HernquistLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	double x = thetaLength/m_angularRadiusScale;

	if (x == 1)
		return m_Sigma_s;

	double x2 = x*x;
	double denom = (x2-1.0);
	denom *= denom;

	return m_Sigma0 * (-3.0 + (2.0 + x2) * F(x))/denom;
}

} // end namespace

