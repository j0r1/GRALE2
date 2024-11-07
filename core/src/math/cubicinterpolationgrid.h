#pragma once

#include "graleconfig.h"
#include "vector2d.h"
#include <errut/booltype.h>
#include <gsl/gsl_interp2d.h>
#include <vector>

namespace grale
{

class GRALE_IMPORTEXPORT CubicInterpolationGrid
{
public:
	CubicInterpolationGrid(Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY); // values need to be set and init needs to be called!
	~CubicInterpolationGrid() { cleanup(); }

	const std::vector<double> &values() const { return m_values; }
	std::vector<double> &values() { return m_values; }

	// Must be called when values changed!
	errut::bool_t init();

	errut::bool_t getValue(Vector2D<double> pos, double *value) const;
	errut::bool_t getGradient(Vector2D<double> pos, Vector2D<double> *pGrad) const;
	errut::bool_t getSecondDerivatives(Vector2D<double> pos, double &dxx, double &dyy, double &dxy) const;
private:
	void cleanup();
	
	const Vector2Dd m_bottomLeft, m_topRight;
	const int m_numX, m_numY;

	gsl_interp2d *m_pInterp;
	gsl_interp_accel *m_pXAccel, *m_pYAccel;
	std::vector<double> m_x, m_y;
	std::vector<double> m_values;
};

} // end namespace

