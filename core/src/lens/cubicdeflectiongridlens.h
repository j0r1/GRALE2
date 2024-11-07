#pragma once

#include "graleconfig.h"
#include "gravitationallens.h"
#include "deflectiongridlens.h"
#include <vector>

namespace grale
{

class CubicInterpolationGrid;

class GRALE_IMPORTEXPORT CubicDeflectionGridLens : public GravitationalLens
{
public:
	CubicDeflectionGridLens();
	~CubicDeflectionGridLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<CubicDeflectionGridLens>(); }

	bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const ;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);

	std::unique_ptr<CubicInterpolationGrid> m_pAxFunction, m_pAyFunction;
	double m_densFactor;
	double m_x0, m_x1, m_y0, m_y1;
};

} // end namespace

