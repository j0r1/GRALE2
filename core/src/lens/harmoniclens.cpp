#include "graleconfig.h"
#include "harmoniclens.h"
#include "constants.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> HarmonicLensParams::createCopy() const
{
	return std::make_unique<HarmonicLensParams>(m_sigma0, m_k, m_l, m_phiX, m_phiY);
}

bool HarmonicLensParams::write(serut::SerializationInterface &si) const
{
	double dparams[5] = { m_sigma0, m_k, m_l, m_phiX, m_phiY };

	if (!si.writeDoubles(dparams, 5))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool HarmonicLensParams::read(serut::SerializationInterface &si)
{
	double dparams[5];

	if (!si.readDoubles(dparams, 5))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_sigma0 = dparams[0];
	m_k = dparams[1];
	m_l = dparams[2];
	m_phiX = dparams[3];
	m_phiY = dparams[4];

	return true;
}


HarmonicLens::HarmonicLens() : GravitationalLens(Harmonic)
{
}

HarmonicLens::~HarmonicLens()
{
}

bool HarmonicLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const HarmonicLensParams *pParams = dynamic_cast<const HarmonicLensParams*>(pLensParams);
	if (pParams == nullptr)
	{
		setErrorString("Invalid parameters");
		return false;
	}

	m_sigma0 = pParams->getDensityScale();
	m_phiX = pParams->getPhiX();
	m_phiY = pParams->getPhiY();
	m_k = pParams->getK();
	m_l = pParams->getL();

	if (m_k == 0 && m_l == 0)
	{
		setErrorString("Both 'k' and 'l' are zero (use MassSheetLens instead)");
		return false;
	}

	double Dd = getLensDistance();
	m_factor = m_sigma0*Dd*(4.0*CONST_PI*CONST_G/(SPEED_C*SPEED_C))*8.0/(m_k*m_k+m_l*m_l);

	return true;
}

bool HarmonicLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	double ax = m_factor*m_k/2.0*SIN(m_k/2.0*theta.getX() + m_phiX)*COS(m_l/2.0*theta.getY() + m_phiY);
	double ay = m_factor*m_l/2.0*COS(m_k/2.0*theta.getX() + m_phiX)*SIN(m_l/2.0*theta.getY() + m_phiY);
	*pAlpha = { ax, ay };
	return true;
}

bool HarmonicLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	axx =  m_factor*m_k*m_k/4.0*COS(m_k/2.0*theta.getX() + m_phiX)*COS(m_l/2.0*theta.getY() + m_phiY);
	ayy =  m_factor*m_l*m_l/4.0*COS(m_k/2.0*theta.getX() + m_phiX)*COS(m_l/2.0*theta.getY() + m_phiY);
	axy = -m_factor*m_k*m_l/4.0*SIN(m_k/2.0*theta.getX() + m_phiX)*SIN(m_l/2.0*theta.getY() + m_phiY);
	return true;
}

double HarmonicLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	return m_sigma0*COS(m_k/2.0*theta.getX() + m_phiX)*COS(m_l/2.0*theta.getY() + m_phiY);
}

bool HarmonicLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	double p = -m_factor*COS(m_k/2.0*theta.getX() + m_phiX)*COS(m_l/2.0*theta.getY() + m_phiY);
	*pPotentialValue = (D_ds/D_s)*p;
	return true;
}

} // end namespace
