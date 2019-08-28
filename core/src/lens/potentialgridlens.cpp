#include "graleconfig.h"
#include "potentialgridlens.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{

GravitationalLensParams *PotentialGridLensParams::createCopy() const
{
	//return new HarmonicLensParams(m_sigma0, m_k, m_l, m_phiX, m_phiY);
	//TODO
	return nullptr;
}

bool PotentialGridLensParams::write(serut::SerializationInterface &si) const
{
	// TODO
	return true;
}

bool PotentialGridLensParams::read(serut::SerializationInterface &si)
{
	return true;
}


PotentialGridLens::PotentialGridLens() : GravitationalLens(PotentialGrid)
{
}

PotentialGridLens::~PotentialGridLens()
{
}

bool PotentialGridLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const PotentialGridLens *pParams = dynamic_cast<const PotentialGridLens*>(pLensParams);
	if (pParams == nullptr)
	{
		setErrorString("Invalid parameters");
		return false;
	}

	// TODO
	//
	return true;
}

bool PotentialGridLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	// TODO
	return true;
}

bool PotentialGridLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	// TODO
	return true;
}

double PotentialGridLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	// TODO
	return 0;
}

bool PotentialGridLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	// TODO
	return true;
}

} // end namespace
