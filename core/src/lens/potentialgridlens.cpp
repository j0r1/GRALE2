#include "graleconfig.h"
#include "potentialgridlens.h"
#include "constants.h"
#include <limits>

using namespace std;
using namespace errut;

namespace grale
{

PotentialGridLensBase::PotentialGridLensBase(double Dd, Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY)
	: CubicInterpolationGrid(bottomLeft, topRight, numX, numY), m_Dd(Dd)
{
}

bool_t PotentialGridLensBase::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	bool_t r;
	if (!(r = getGradient(theta, pAlpha)))
		return "Error getting derivative of potential: " + r.getErrorString();
	return true;
}

bool_t PotentialGridLensBase::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	bool_t r;

	if (!(r = getSecondDerivatives(theta, axx, ayy, axy)))
			return "Error getting second derivative of potential: " + r.getErrorString();

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
	bool_t r;
	if (!(r = getValue(theta, pPotentialValue)))
		return "Unable to get projected potential at requested point: " + r.getErrorString();

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
	if (numX*numY != (int)values.size())
	{
		setErrorString("Number of potential values doesn't match grid dimensions");
		return false;
	}

	if (numX < (int)gsl_interp2d_type_min_size(gsl_interp2d_bicubic) ||
		numY < (int)gsl_interp2d_type_min_size(gsl_interp2d_bicubic))
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
