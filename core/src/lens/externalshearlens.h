#pragma once

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

class GRALE_IMPORTEXPORT ExternalShearLensParams : public GravitationalLensParams
{
public:
	ExternalShearLensParams()											{ }
	ExternalShearLensParams(double shearSize, double shearAngleDeg)		{ m_shearSize = shearSize; m_shearAngle = shearAngleDeg; }

	double getShearSize() const											{ return m_shearSize; }
	double getShearAngleDegrees() const									{ return m_shearAngle; }
	
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double m_shearSize = 0;
	double m_shearAngle = 0;
};

class GRALE_IMPORTEXPORT ExternalShearLens : public GravitationalLens
{
public:
	ExternalShearLens();
	~ExternalShearLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<ExternalShearLens>(); }
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
	std::vector<CLFloatParamInfo> getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const override;
	std::unique_ptr<GravitationalLensParams> createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const override;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double m_shearSize = 0;
	double m_shearAngleRad = 0;
	double m_gamma1 = 0, m_gamma2 = 0;
};

} // end namespace
