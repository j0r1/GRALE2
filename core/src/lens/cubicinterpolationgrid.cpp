#include "graleconfig.h"
#include "cubicinterpolationgrid.h"
#include <limits>

using namespace std;
using namespace errut;

namespace grale
{

CubicInterpolationGrid::CubicInterpolationGrid(Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY)
	: m_bottomLeft(bottomLeft), m_topRight(topRight), m_numX(numX), m_numY(numY)
{
	m_values.resize(numX*numY);
	m_pInterp = nullptr;
	m_pXAccel = nullptr;
	m_pYAccel = nullptr;

	for (int x = 0 ; x < numX ; x++)
		m_x.push_back((topRight.getX()-bottomLeft.getX())/(double)(numX-1) * x + bottomLeft.getX());

	for (int y = 0 ; y < numY ; y++)
		m_y.push_back((topRight.getY()-bottomLeft.getY())/(double)(numY-1) * y + bottomLeft.getY());
}

errut::bool_t CubicInterpolationGrid::init()
{
	cleanup();

	if (m_bottomLeft.getX() >= m_topRight.getX() || m_bottomLeft.getY() >= m_topRight.getY())
		return "Corners need to be ordered correctly";

	if (m_numX*m_numY != (int)m_values.size())
		return "Number of potential values doesn't match grid dimensions";

	if (m_numX < (int)gsl_interp2d_type_min_size(gsl_interp2d_bicubic) ||
		m_numY < (int)gsl_interp2d_type_min_size(gsl_interp2d_bicubic))
		return "Not enough points in grid in at least one dimension";

	m_pInterp = gsl_interp2d_alloc(gsl_interp2d_bicubic, m_numX, m_numY);
	if (!m_pInterp)
		return "Unable to allocate GSL interpolation workspace";

	m_pXAccel = gsl_interp_accel_alloc();
	m_pYAccel = gsl_interp_accel_alloc();

	int err = gsl_interp2d_init(m_pInterp, m_x.data(), m_y.data(), m_values.data(), m_numX, m_numY);
	if (err < 0)
		return "Unable to initialize GSL interpolation workspace: error " + to_string(err);

	return true;
}

void CubicInterpolationGrid::cleanup()
{
	if (m_pXAccel)
		gsl_interp_accel_free(m_pXAccel);
	if (m_pYAccel)
		gsl_interp_accel_free(m_pYAccel);
	if (m_pInterp)
		gsl_interp2d_free(m_pInterp);

	m_pInterp = nullptr;
	m_pXAccel = nullptr;
	m_pYAccel = nullptr;	
}

bool_t CubicInterpolationGrid::getGradient(Vector2D<double> pos, Vector2D<double> *pGrad) const
{
	double dx = 0, dy = 0;

	auto getDeriv = [this, pos](
			int (*functionName)(const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d),
			double &result) -> bool_t
	{
		int err = functionName(m_pInterp, m_x.data(), m_y.data(), m_values.data(),
				            pos.getX(), pos.getY(), m_pXAccel, m_pYAccel, &result);
		if (err < 0)
			return "Error getting derivative of interpolated grid: GSL error code " + to_string(err);
		return true;
	};

	bool_t r = getDeriv(gsl_interp2d_eval_deriv_x_e, dx);
	if (!r)
		return r;
	r = getDeriv(gsl_interp2d_eval_deriv_y_e, dy);
	if (!r)
		return r;

	*pGrad = Vector2Dd(dx, dy);
	return true;
}

bool_t CubicInterpolationGrid::getSecondDerivatives(Vector2D<double> pos, double &dxx, double &dyy, double &dxy) const
{
	auto getDerivDeriv = [this, pos](
			int (*functionName)(const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d),
			double &result) -> bool_t
	{
		int err = functionName(m_pInterp, m_x.data(), m_y.data(), m_values.data(), 
				               pos.getX(), pos.getY(), m_pXAccel, m_pYAccel, &result);
		if (err < 0)
			return "Error getting second derivative of interpolated grid: GSL error code " + to_string(err);
		return true;
	};

	bool_t r = getDerivDeriv(gsl_interp2d_eval_deriv_xx_e, dxx);
	if (!r)
		return r;

	r = getDerivDeriv(gsl_interp2d_eval_deriv_yy_e, dyy);
	if (!r)
		return r;

	r = getDerivDeriv(gsl_interp2d_eval_deriv_xy_e, dxy);
	if (!r)
		return r;

	return true;
}

bool_t CubicInterpolationGrid::getValue(Vector2D<double> pos, double *pValue) const
{
	double result;
	int err = gsl_interp2d_eval_e(m_pInterp, m_x.data(), m_y.data(), m_values.data(),
	                              pos.getX(), pos.getY(), m_pXAccel, m_pYAccel, &result);
	if (err < 0)
		return "Unable to get interpolated value at requested point: GSL error code " + to_string(err);

	*pValue = result;
	return true;
}

} // end namespace
