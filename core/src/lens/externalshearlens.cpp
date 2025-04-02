#include "graleconfig.h"
#include "externalshearlens.h"
#include "constants.h"

using namespace std;

namespace grale
{

std::unique_ptr<GravitationalLensParams> ExternalShearLensParams::createCopy() const
{
	return make_unique<ExternalShearLensParams>(m_shearSize, m_shearAngle);
}

bool ExternalShearLensParams::write(serut::SerializationInterface &si) const
{
	if (!si.writeDouble(m_shearSize) || !si.writeDouble(m_shearAngle))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool ExternalShearLensParams::read(serut::SerializationInterface &si)
{
	if (!si.readDouble(&m_shearSize) || !si.readDouble(&m_shearAngle))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

ExternalShearLens::ExternalShearLens() : GravitationalLens(ExternalShear)
{
}

ExternalShearLens::~ExternalShearLens()
{
}

bool ExternalShearLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const ExternalShearLensParams *pParams = dynamic_cast<const ExternalShearLensParams *>(pLensParams);
	if (!pParams)
	{
		setErrorString("Parameters are not of type 'ExternalShearLensParams'");
		return false;
	}

	m_shearSize = pParams->getShearSize();
	m_shearAngleRad = (pParams->getShearAngleDegrees()/180.0)*CONST_PI;

	m_gamma1 = m_shearSize * std::cos(2.0*m_shearAngleRad);
	m_gamma2 = m_shearSize * std::sin(2.0*m_shearAngleRad);

	return true;
}

bool ExternalShearLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	double ax = m_gamma1*theta.getX() + m_gamma2*theta.getY();
	double ay = m_gamma2*theta.getX() - m_gamma1*theta.getY();
	*pAlpha = Vector2Dd(ax, ay);
	return true;
}

bool ExternalShearLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	axx = m_gamma1;
	ayy = -m_gamma1;
	axy = m_gamma2;
	return true;
}

double ExternalShearLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	return 0;
}

bool ExternalShearLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	*pPotentialValue = 0.5*m_gamma1*(theta.getX()*theta.getX() - theta.getY()*theta.getY()) + m_gamma2*theta.getX()*theta.getY();
	return true;
}

bool ExternalShearLens::getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const
{
	// TODO: there's no real scale here, what to use?
	*pDeflectionScale = ANGLE_ARCSEC;
	*pPotentialScale = ANGLE_ARCSEC*ANGLE_ARCSEC;
	return true;
}

bool ExternalShearLens::getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const
{
	*pNumIntParams = 0;
	*pNumFloatParams = 2; // size and angle
	return true;
}

bool ExternalShearLens::getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const
{
	pFloatParams[0] = (float)m_shearSize;
	pFloatParams[1] = (float)m_shearAngleRad;
	return true;
}

std::string ExternalShearLens::getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives, bool potential) const
{
	string program;
	subRoutineName = "clExternalShearProgram";

	program += R"XYZ(
LensQuantities clExternalShearProgram(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
{
	float shearSize = pFloatParams[0];
	float shearAngle = pFloatParams[1];
	float gamma1 = shearSize*cos(2.0*shearAngle);
	float gamma2 = shearSize*sin(2.0*shearAngle);

	LensQuantities r;
	r.alphaX = gamma1*coord.x + gamma2*coord.y;
	r.alphaY = gamma2*coord.x - gamma1*coord.y;
)XYZ";
	if (potential)
		program += R"XYZ(
	r.potential = 0.5*gamma1*(coord.x*coord.x - coord.y*coord.y) + gamma2*coord.x*coord.y;
)XYZ";
	if (derivatives)
		program += R"XYZ(
	r.axx = gamma1;
	r.ayy = -gamma1;
	r.axy = gamma2;
)XYZ";
	program += R"XYZ(
	return r;
}
)XYZ";
	return program;
}

std::vector<CLFloatParamInfo> ExternalShearLens::getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const
{
	return {
		{ .name = "shearsize", .offset = 0, .scaleFactor = 1.0, .hardMin = 0 },
		{ .name = "shearangle", .offset = 1, .scaleFactor = 180.0/CONST_PI },
	};
}

std::unique_ptr<GravitationalLensParams> ExternalShearLens::createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const
{
	double shearSize = (double)pFloatParams[0];
	double shearAngleDeg = (((double)pFloatParams[1])/CONST_PI)*180.0;

	return make_unique<ExternalShearLensParams>(shearSize, shearAngleDeg);
}


} // end namespace
