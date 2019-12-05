#include "graleconfig.h"
#include "circularpieceslens.h"
#include "constants.h"
#include <limits>

#include "debugnew.h"

using namespace std;

namespace grale
{

CircularPiecesLensParams::CircularPiecesLensParams(const vector<CircularPieceInfo> &pieces, const vector<double> &coeffs)
	: m_pieces(pieces),
	  m_coeffs(coeffs)
{
}

GravitationalLensParams *CircularPiecesLensParams::createCopy() const
{
	vector<CircularPieceInfo> pieces;
	for (auto &p : m_pieces)
	{
		// Create actual copy of lens model
		auto lens = p.getLens();
		if (lens.get() == nullptr)
		{
			setErrorString("A piece does not contain a lens model");
			return nullptr;
		}

		shared_ptr<GravitationalLens> newLens(lens->createCopy());
		if (newLens.get() == nullptr)
		{
			setErrorString("Unable to create copy of a lens: " + lens->getErrorString());
			return nullptr;
		}

		pieces.push_back({ newLens, p.getStartRadius(), p.getEndRadius(), 
				           p.getPotentialScale(), p.getPotentialOffset() });
	}

	CircularPiecesLensParams *pNewParams = new CircularPiecesLensParams();
	std::swap(pieces, pNewParams->m_pieces);
	pNewParams->m_coeffs = m_coeffs;
	return pNewParams;
}

bool CircularPiecesLensParams::write(serut::SerializationInterface &si) const
{
	int32_t numPieces = m_pieces.size();
	if (!si.writeInt32(numPieces))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	for (auto &p : m_pieces)
	{
		double v[4];

		v[0] = p.getStartRadius();
		v[1] = p.getEndRadius();
		v[2] = p.getPotentialScale();
		v[3] = p.getPotentialOffset();

		auto lens = p.getLens();
		if (lens.get() == nullptr)
		{
			setErrorString("A piece does not contain a lens model");
			return false;
		}

		if (!lens->write(si) || !si.writeDoubles(v, 4))
		{
			setErrorString(si.getErrorString());
			return false;
		}
	}

	int32_t numCoeffs = m_coeffs.size();
	if (!si.writeInt32(numCoeffs) || !si.writeDoubles(m_coeffs))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	return true;
}

bool CircularPiecesLensParams::read(serut::SerializationInterface &si)
{
	int32_t numPieces;
	if (!si.readInt32(&numPieces))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_pieces.clear();

	for (int32_t i = 0 ; i < numPieces ; i++)
	{
		GravitationalLens *pLens = nullptr;
		string errStr;

		if (!GravitationalLens::read(si, &pLens, errStr))
		{
			setErrorString("Couldn't read a lens model");
			return false;
		}

		shared_ptr<GravitationalLens> lens(pLens);
		double v[4];

		if (!si.readDoubles(v, 4))
		{
			setErrorString(si.getErrorString());
			return false;
		}

		m_pieces.push_back({ lens, v[0], v[1], v[2], v[3] });
	}

	int32_t numCoeffs;
	if (!si.readInt32(&numCoeffs))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_coeffs.resize(numCoeffs);
	if (!si.readDoubles(m_coeffs))
	{
		setErrorString(si.getErrorString());
		return false;
	}

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
	const CircularPiecesLensParams *pParams = dynamic_cast<const CircularPiecesLensParams*>(pLensParams);
	if (pParams == nullptr)
	{
		setErrorString("Invalid parameters");
		return false;
	}

	size_t numPieces = pParams->getPiecesInfo().size();
	if (numPieces == 0)
	{
		setErrorString("No pieces were specified");
		return false;
	}

	// Create copy to get own instances of the lenses
	CircularPiecesLensParams *pCopy = static_cast<CircularPiecesLensParams*>(pParams->createCopy());
	shared_ptr<CircularPiecesLensParams> params(pCopy);
	auto pieces = params->getPiecesInfo();

	if (pieces[0].getStartRadius() != 0)
	{
		setErrorString("Start radius of first piece must be zero");
		return false;
	}

	if (pieces[numPieces-1].getEndRadius() != numeric_limits<double>::infinity())
	{
		setErrorString("End radius of last piece must be infinity");
		return false;
	}

	for (size_t i = 0 ; i < pieces.size() ; i++)
	{
		if (pieces[i].getStartRadius() > pieces[i].getEndRadius())
		{
			setErrorString("Start radius is greater than end radius for piece " + to_string(i+1));
			return false;
		}

		if (i+1 < pieces.size())
		{
			if (pieces[i].getEndRadius() > pieces[i+1].getStartRadius())
			{
				setErrorString("End radius for piece " + to_string(i+1) + " is larger than start radius of next piece");
				return false;
			}
		}

		auto &p = pieces[i];
		m_lenses.push_back(p.getLens());
		m_startRadius.push_back(p.getStartRadius());
		m_endRadius.push_back(p.getEndRadius());
		m_scale.push_back(p.getPotentialScale());
		m_potentialOffset.push_back(p.getPotentialOffset());
	}

	auto getDerivsCoeffs = [](const vector<double> &coeffs) -> vector<double>
	{
		vector<double> result;
		for (size_t i = 1 ; i < coeffs.size() ; i++)
			result.push_back((double)i * coeffs[i]);
		return result;
	};

	m_fCoeffs = pParams->getInterpolationFunctionCoefficients();
	m_fAccCoeffs = getDerivsCoeffs(m_fCoeffs);
	m_fAccAccCoeffs = getDerivsCoeffs(m_fAccCoeffs);

	return true;
}

void CircularPiecesLens::getLenses(Vector2D<double> theta, int &idx1, int &idx2, double &x, double &diffTheta) const
{
	double r = theta.getLength();
	size_t numPieces = m_lenses.size();

	idx1 = -1;
	idx2 = -1;
	x = 0;
	diffTheta = 0;

	if (r >= m_startRadius[numPieces-1]) // Must be last lens only
	{
		idx1 = (int)numPieces-1;
		return;
	}

	for (int i = 0 ; i < (int)m_lenses.size() - 1 ; i++)
	{
		if (r <= m_endRadius[i]) // Lies completely within this lens
		{
			idx1 = i;
			return;
		}

		if (r <= m_startRadius[i+1]) // Lies in between i and i+1
		{
			idx1 = i;
			idx2 = i+1;
			diffTheta = m_startRadius[i+1]-m_endRadius[i];
			x = (r-m_endRadius[i])/diffTheta;
			return;
		}
	}

	setErrorString("Internal error: should not happen");
}

//#define NUMERICDIFFERENTIATION
#ifdef NUMERICDIFFERENTIATION
static double D = 0.0001*ANGLE_ARCSEC; // TODO: remove this again after testing
#endif // NUMERICDIFFERENTIATION

bool CircularPiecesLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
#ifdef NUMERICDIFFERENTIATION
	double p0, px, py;
	getProjectedPotential(1, 1, theta, &p0);
	getProjectedPotential(1, 1, theta + Vector2Dd(D, 0), &px);
	getProjectedPotential(1, 1, theta + Vector2Dd(0, D), &py);

	*pAlpha = Vector2Dd((px-p0)/D, (py-p0)/D);
#else
	int idx1, idx2;
	double x, diffTheta;

	getLenses(theta, idx1, idx2, x, diffTheta);
	if (idx1 < 0)
		return false;

	auto getSub = [this,theta](int idx, Vector2Dd &alpha, double &v)
	{
		if (!m_lenses[idx]->getAlphaVector(theta, &alpha))
		{
			setErrorString("Couldn't get piese lens deflection angle: " + m_lenses[idx]->getErrorString());
			return false;
		}
		alpha *= m_scale[idx];

		if (!m_lenses[idx]->getProjectedPotential(1, 1, theta, &v))
		{
			setErrorString("Couldn't get piece lens potential: " + m_lenses[idx]->getErrorString());
			return false;
		}

		v += m_potentialOffset[idx];
		v *= m_scale[idx];
		return true;
	};

	Vector2Dd alpha1, alpha2;
	double v1, v2;
	if (!getSub(idx1, alpha1, v1))
		return false;

	if (idx2 < 0) // Just one lens, scaled potential
	{
		*pAlpha = alpha1;
		return true;
	}

	// Need to interpolate
	if (!getSub(idx2, alpha2, v2))
		return false;

	double f = getF(x);
	double facc = getFDeriv(x);

	double thetaLength = theta.getLength();
	Vector2Dd result = theta/thetaLength;

	result *= (facc/diffTheta)*(v1-v2);
	result += alpha1*f;
	result += alpha2*(1.0-f);
	*pAlpha = result;
#endif
	return true;
}

bool CircularPiecesLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
#ifdef NUMERICDIFFERENTIATION
	Vector2Dd a0, ax, ay;

	getAlphaVector(theta, &a0);
	getAlphaVector(theta + Vector2Dd(D, 0), &ax);
	getAlphaVector(theta + Vector2Dd(0, D), &ay);

	axx = (ax.getX()-a0.getX())/D;
	ayy = (ay.getY()-a0.getY())/D;
	axy = 0.5*(ax.getY() - a0.getY())/D + 0.5*(ay.getX() - a0.getX())/D;
#else
	int idx1, idx2;
	double x, diffTheta;

	getLenses(theta, idx1, idx2, x, diffTheta);
	if (idx1 < 0)
		return false;

	auto getSub = [this,theta](int idx, double &axx, double &ayy, double &axy, Vector2Dd &alpha, double &v)
	{
		if (!m_lenses[idx]->getAlphaVectorDerivatives(theta, axx, ayy, axy))
		{
			setErrorString("Coudln't get piece lens deflection angle derivatives: " + m_lenses[idx]->getErrorString());
			return false;
		}
		axx *= m_scale[idx];
		ayy *= m_scale[idx];
		axy *= m_scale[idx];

		if (!m_lenses[idx]->getAlphaVector(theta, &alpha))
		{
			setErrorString("Couldn't get piece lens deflection angle: " + m_lenses[idx]->getErrorString());
			return false;
		}
		alpha *= m_scale[idx];

		if (!m_lenses[idx]->getProjectedPotential(1, 1, theta, &v))
		{
			setErrorString("Couldn't get piece lens potential: " + m_lenses[idx]->getErrorString());
			return false;
		}

		v += m_potentialOffset[idx];
		v *= m_scale[idx];
		return true;
	};

	double axx1, ayy1, axy1, axx2, ayy2, axy2;
	Vector2Dd alpha1, alpha2;
	double v1, v2;

	if (!getSub(idx1, axx1, ayy1, axy1, alpha1, v1))
		return false;

	if (idx2 < 0) // Just one lens
	{
		axx = axx1;
		ayy = ayy1;
		axy = axy1;
		return true;
	}

	// Interpolation
	if (!getSub(idx2, axx2, ayy2, axy2, alpha2, v2))
		return false;

	double f = getF(x);
	double facc = getFDeriv(x);
	double faccacc = getFDerivDeriv(x);
	double thetaLength = theta.getLength();
	Vector2Dd unit = theta/thetaLength;
	Vector2Dd alphaDiff = alpha1 - alpha2;

	double faccScaled = facc/diffTheta;
	double faccaccScaled2 = faccacc/(diffTheta*diffTheta);

	axx = axx1*f + axx2*(1.0-f) + 2.0*faccScaled*unit.getX()*(alpha1.getX() - alpha2.getX())
		+ (v1-v2)*(faccaccScaled2*unit.getX()*unit.getX() + (faccScaled/thetaLength)*unit.getY()*unit.getY());

	ayy = ayy1*f + ayy2*(1.0-f) + 2.0*faccScaled*unit.getY()*(alpha1.getY() - alpha2.getY())
		+ (v1-v2)*(faccaccScaled2*unit.getY()*unit.getY() + (faccScaled/thetaLength)*unit.getX()*unit.getX());

	axy = axy1*f + axy2*(1.0-f) + alphaDiff.getX()*faccScaled*unit.getY() 
	                            + alphaDiff.getY()*faccScaled*unit.getX()
	    + (v1-v2)*(faccaccScaled2 - faccScaled/thetaLength)*unit.getX()*unit.getY();
#endif // NUMERICDIFFERENTIATION
	return true;
}

double CircularPiecesLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double axx, ayy, axy;

	getAlphaVectorDerivatives(theta, axx, ayy, axy);
	double kappa = 0.5*(axx + ayy);
	double sigmaCrit = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*getLensDistance());
	return kappa*sigmaCrit;
}

bool CircularPiecesLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	int idx1, idx2;
	double x, diffTheta;

	getLenses(theta, idx1, idx2, x, diffTheta);
	if (idx1 < 0)
		return false;

	auto getSub = [this,D_s,D_ds,theta](int idx, double &v)
	{
		if (!m_lenses[idx]->getProjectedPotential(D_s, D_ds, theta, &v))
		{
			setErrorString("Couldn't get piece lens potential: " + m_lenses[idx]->getErrorString());
			return false;
		}

		v += m_potentialOffset[idx];
		v *= m_scale[idx];
		return true;
	};

	double v1;
	if (!getSub(idx1, v1))
		return false;

	if (idx2 < 0) // Just a single lens
	{
		*pPotentialValue = v1;
		return true;
	}

	double v2;
	if (!getSub(idx2, v2))
		return false;

	double f = getF(x);
	*pPotentialValue = f*v1 + (1.0-f)*v2;
	return true;
}

static inline double calcFunction(const vector<double> &coeffs, double x)
{
	assert(x >= 0 && x <= 1);
	double factor = 1.0;
	double result = 0;
	for (auto a : coeffs)
	{
		result += a*factor;
		factor *= x;
	}
	return result;
}

inline double CircularPiecesLens::getF(double x) const
{
	return calcFunction(m_fCoeffs, x);
}

inline double CircularPiecesLens::getFDeriv(double x) const
{
	return calcFunction(m_fAccCoeffs, x);
}

inline double CircularPiecesLens::getFDerivDeriv(double x) const
{
	return calcFunction(m_fAccAccCoeffs, x);
}


} // end namespace
