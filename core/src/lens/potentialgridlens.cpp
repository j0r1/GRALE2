#include "graleconfig.h"
#include "potentialgridlens.h"
#include "constants.h"

using namespace std;
using namespace errut;

namespace grale
{

PotentialGridLensBase::PotentialGridLensBase(double Dd, Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY)
	: m_Dd(Dd), m_bottomLeft(bottomLeft), m_topRight(topRight), m_numX(numX), m_numY(numY)
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

errut::bool_t PotentialGridLensBase::init()
{
	cleanup();

	if (m_bottomLeft.getX() >= m_topRight.getX() || m_bottomLeft.getY() >= m_topRight.getY())
		return "Corners need to be ordered correctly";

	if (m_numX*m_numY != m_values.size())
		return "Number of potential values doesn't match grid dimensions";

	if (m_numX < gsl_interp2d_type_min_size(gsl_interp2d_bicubic) ||
		m_numY < gsl_interp2d_type_min_size(gsl_interp2d_bicubic))
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

void PotentialGridLensBase::cleanup()
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

bool_t PotentialGridLensBase::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	double ax = 0, ay = 0;
	auto getDeriv = [this, theta](
			int (*functionName)(const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d),
			double &result) -> bool_t
	{
		int err = functionName(m_pInterp, m_x.data(), m_y.data(), m_values.data(),
				            theta.getX(), theta.getY(), m_pXAccel, m_pYAccel, &result);
		if (err < 0)
			return "Error getting derivative of potential: GSL error code " + to_string(err);
		return true;
	};

	bool_t r = getDeriv(gsl_interp2d_eval_deriv_x_e, ax);
	if (!r)
		return r;
	r = getDeriv(gsl_interp2d_eval_deriv_y_e, ay);
	if (!r)
		return r;

	*pAlpha = Vector2Dd(ax, ay);
	return true;
}

bool_t PotentialGridLensBase::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	auto getDerivDeriv = [this, theta](
			int (*functionName)(const gsl_interp2d * interp, const double xa[], const double ya[], const double za[], const double x, const double y, gsl_interp_accel * xacc, gsl_interp_accel * yacc, double * d),
			double &result) -> bool_t
	{
		int err = functionName(m_pInterp, m_x.data(), m_y.data(), m_values.data(), 
				               theta.getX(), theta.getY(), m_pXAccel, m_pYAccel, &result);
		if (err < 0)
			return "Error getting second derivative of potential: GSL error code " + to_string(err);
		return true;
	};

	bool_t r = getDerivDeriv(gsl_interp2d_eval_deriv_xx_e, axx);
	if (!r)
		return r;

	r = getDerivDeriv(gsl_interp2d_eval_deriv_yy_e, ayy);
	if (!r)
		return r;

	r = getDerivDeriv(gsl_interp2d_eval_deriv_xy_e, axy);
	if (!r)
		return r;

	return true;
}

bool_t PotentialGridLensBase::getSurfaceMassDensity(Vector2D<double> theta, double &dens) const
{
	double axx, ayy, axy;

	bool_t r = getAlphaVectorDerivatives(theta, axx, ayy, axy);
	if (!r)
		return r;

	double kappa = 0.5*(axx + ayy);
	double sigmaCrit = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*m_Dd);
	dens = kappa*sigmaCrit;
	return true;
}

bool_t PotentialGridLensBase::getProjectedPotential(Vector2D<double> theta, double *pPotentialValue) const
{
	double result;
	int err = gsl_interp2d_eval_e(m_pInterp, m_x.data(), m_y.data(), m_values.data(),
	                              theta.getX(), theta.getY(), m_pXAccel, m_pYAccel, &result);
	if (err < 0)
		return "Unable to get projected potential at requested point: GSL error code " + to_string(err);

	*pPotentialValue = result;
	return true;
}

PotentialGridLensParams::PotentialGridLensParams()
{
}

PotentialGridLensParams::PotentialGridLensParams(Vector2Dd bottomLeft, Vector2Dd topRight, 
                                                 const std::vector<double> &values, int numX, int numY)
	: m_bottomLeft(bottomLeft),
	  m_topRight(topRight),
	  m_values(values),
	  m_numX(numX),
	  m_numY(numY)
{
}

PotentialGridLensParams::~PotentialGridLensParams()
{
}

std::unique_ptr<GravitationalLensParams> PotentialGridLensParams::createCopy() const
{
	return std::make_unique<PotentialGridLensParams>(m_bottomLeft, m_topRight, m_values, m_numX, m_numY);
}

bool PotentialGridLensParams::write(serut::SerializationInterface &si) const
{
	double region[4] = { m_bottomLeft.getX(), m_bottomLeft.getY(), m_topRight.getX(), m_topRight.getY() };
	int32_t numNXY[3] = { (int32_t)m_values.size(), m_numX, m_numY };

	if (!si.writeDoubles(region, 4) ||
		!si.writeInt32s(numNXY, 3) ||
		!si.writeDoubles(m_values))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool PotentialGridLensParams::read(serut::SerializationInterface &si)
{
	double region[4];
	int32_t nxy[3];

	if (!si.readDoubles(region, 4) || !si.readInt32s(nxy, 3))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_bottomLeft = Vector2Dd(region[0], region[1]);
	m_topRight = Vector2Dd(region[2], region[3]);
	m_values.resize(nxy[0]);
	m_numX = nxy[1];
	m_numY = nxy[2];

	if (!si.readDoubles(m_values))
	{
		setErrorString(si.getErrorString());
		return false;
	}
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
	const PotentialGridLensParams *pParams = dynamic_cast<const PotentialGridLensParams*>(pLensParams);
	if (pParams == nullptr)
	{
		setErrorString("Invalid parameters");
		return false;
	}

	Vector2Dd bottomLeft = pParams->getBottomLeft();
	Vector2Dd topRight = pParams->getTopRight();

	if (bottomLeft.getX() >= topRight.getX() || bottomLeft.getY() >= topRight.getY())
	{
		setErrorString("Corners need to be ordered correctly");
		return false;
	}

	const vector<double> &values = pParams->getValues();
	int numX = pParams->getNumX();
	int numY = pParams->getNumY();
	if (numX*numY != values.size())
	{
		setErrorString("Number of potential values doesn't match grid dimensions");
		return false;
	}

	if (numX < gsl_interp2d_type_min_size(gsl_interp2d_bicubic) ||
		numY < gsl_interp2d_type_min_size(gsl_interp2d_bicubic))
	{
		setErrorString("Not enough points in grid in at least one dimension");
		return false;
	}

	m_pBase = make_unique<PotentialGridLensBase>(getLensDistance(), bottomLeft, topRight, numX, numY);
	m_pBase->values() = values;

	bool_t r = m_pBase->init();
	if (!r)
	{
		m_pBase = nullptr;
		setErrorString(r.getErrorString());
		return false;
	}

	return true;
}

bool PotentialGridLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	bool_t r = m_pBase->getAlphaVector(theta, pAlpha);
	if (!r)
	{
		setErrorString(r.getErrorString());
		return false;
	}
	return true;
}

bool PotentialGridLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	bool_t r = m_pBase->getAlphaVectorDerivatives(theta, axx, ayy, axy);
	if (!r)
	{
		setErrorString(r.getErrorString());
		return false;
	}
	return true;
}

// TODO: make this the default in GravitationalLens base class?
double PotentialGridLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double dens;
	bool_t r = m_pBase->getSurfaceMassDensity(theta, dens);
	if (!r)
	{
		setErrorString(r.getErrorString());
		return numeric_limits<double>::quiet_NaN();
	}
	return dens;
}

bool PotentialGridLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	double v;
	bool_t r = m_pBase->getProjectedPotential(theta, &v);
	if (!r)
	{
		setErrorString(r.getErrorString());
		return false;
	}

	*pPotentialValue = v*(D_ds/D_s);
	return true;
}

} // end namespace
