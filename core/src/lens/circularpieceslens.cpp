#include "graleconfig.h"
#include "circularpieceslens.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{

GravitationalLensParams *CircularPiecesLensParams::createCopy() const
{
	//return new HarmonicLensParams(m_sigma0, m_k, m_l, m_phiX, m_phiY);
	//TODO
	return nullptr;
}

bool CircularPiecesLensParams::write(serut::SerializationInterface &si) const
{
	// TODO
	return true;
}

bool CircularPiecesLensParams::read(serut::SerializationInterface &si)
{
	return true;
}


CircularPiecesLens::CircularPiecesLens() : GravitationalLens(CircularPieces)
{
}

CircularPiecesLens::~CircularPiecesLens()
{
}

bool CircularPiecesLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const CircularPiecesLens *pParams = dynamic_cast<const CircularPiecesLens*>(pLensParams);
	if (pParams == nullptr)
	{
		setErrorString("Invalid parameters");
		return false;
	}

	// TODO
	//
	return true;
}

bool CircularPiecesLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	// TODO
	return true;
}

bool CircularPiecesLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	// TODO
	return true;
}

double CircularPiecesLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	// TODO
	return 0;
}

bool CircularPiecesLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	// TODO
	return true;
}

} // end namespace
