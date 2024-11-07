#include "graleconfig.h"
#include "cubicdeflectiongridlens.h"
#include "cubicinterpolationgrid.h"
#include "constants.h"
#include <cmath>
#include <iostream>
#include <limits>

using namespace std;
using namespace errut;

namespace grale
{

// Same parameters as DeflectionGridLens, different lens

CubicDeflectionGridLens::CubicDeflectionGridLens() : GravitationalLens(CubicDeflectionGrid)
{
}

CubicDeflectionGridLens::~CubicDeflectionGridLens()
{
}

bool CubicDeflectionGridLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const DeflectionGridLensParams *pParams = dynamic_cast<const DeflectionGridLensParams *>(pLensParams);

	if (pParams == 0)
	{
		setErrorString("Parameters are not of type 'DeflectionGridLensParams'");
		return false;
	}


	Vector2D<double> bottomLeft = pParams->getBottomLeft();
	Vector2D<double> topRight = pParams->getTopRight();

	m_pAxFunction = make_unique<CubicInterpolationGrid>(bottomLeft, topRight, pParams->getWidth(), pParams->getHeight());
	m_pAyFunction = make_unique<CubicInterpolationGrid>(bottomLeft, topRight, pParams->getWidth(), pParams->getHeight());
	m_pAxFunction->values() = pParams->getAlphaX();
	m_pAyFunction->values() = pParams->getAlphaY();

	bool_t r;
	if (!(r = m_pAyFunction->init()) || !(r = m_pAxFunction->init()))
	{
		setErrorString(r.getErrorString());
		return false;
	}

	m_densFactor = (SPEED_C*SPEED_C)/(getLensDistance()*8.0*CONST_PI*CONST_G);

	m_x0 = bottomLeft.getX();
	m_x1 = topRight.getX();
	m_y0 = bottomLeft.getY();
	m_y1 = topRight.getY();

	if (m_x0 > m_x1) std::swap(m_x0, m_x1);
	if (m_y0 > m_y1) std::swap(m_y0, m_y1);

	return true;
}

bool CubicDeflectionGridLens::getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const
{
	if (theta.getX() < m_x0 || theta.getX() > m_x1 ||
		theta.getY() < m_y0 || theta.getY() > m_y1)
	{
		setErrorString("Theta position doesn't lie inside area of deflection grid");
		return false;
	}

	bool_t r;
	double alphaX, alphaY;
	if (!(r = m_pAxFunction->getValue(theta, &alphaX)) || !(r = m_pAyFunction->getValue(theta, &alphaY)))
	{
		setErrorString(r.getErrorString());
		return false;
	}

	*pAlpha = Vector2D<double>(alphaX, alphaY);
	return true;
}

double CubicDeflectionGridLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	if (theta.getX() < m_x0 || theta.getX() > m_x1 ||
		theta.getY() < m_y0 || theta.getY() > m_y1)
		return std::numeric_limits<double>::quiet_NaN();

	Vector2Dd gradAx, gradAy;
	bool_t r;
	if (!(r = m_pAxFunction->getGradient(theta, &gradAx)) || !(r = m_pAyFunction->getGradient(theta, &gradAy)))
		return std::numeric_limits<double>::quiet_NaN();

	double daxx = gradAx.getX();
	double dayy = gradAy.getY();

	return m_densFactor*(daxx + dayy);
}

bool CubicDeflectionGridLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	if (theta.getX() < m_x0 || theta.getX() > m_x1 ||
		theta.getY() < m_y0 || theta.getY() > m_y1)
	{
		setErrorString("Theta position doesn't lie inside area of deflection grid");
		return false;
	}

	Vector2Dd gradAx, gradAy;
	bool_t r;
	if (!(r = m_pAxFunction->getGradient(theta, &gradAx)) || !(r = m_pAyFunction->getGradient(theta, &gradAy)))
	{
		setErrorString(r.getErrorString());
		return false;
	}

	axx = gradAx.getX();
	ayy = gradAy.getY();
	axy = 0.5*(gradAx.getY() + gradAy.getX());

	return true;
}

bool CubicDeflectionGridLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	setErrorString("getProjectedPotential is not supported for CubicDeflectionGridLens");
	return false;
}

} // end namespace

