#pragma once

#include "graleconfig.h"
#include "gravitationallens.h"
#include <vector>
#include <memory>
#include <assert.h>

namespace grale
{

class GRALE_IMPORTEXPORT MultiPlaneContainerParams : public GravitationalLensParams
{
public:
	MultiPlaneContainerParams();
	~MultiPlaneContainerParams();

	bool add(std::shared_ptr<GravitationalLens> lens, double z);
	bool add(GravitationalLens *pLens, double z);

	int getNumberOfLenses() const { return m_lensesAndRedshifts.size(); }
	const GravitationalLens &getLens(int idx) const;
	double getRedshift(int idx) const;

	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	std::vector<std::pair<std::shared_ptr<GravitationalLens>, double>> m_lensesAndRedshifts;
};

inline const GravitationalLens &MultiPlaneContainerParams::getLens(int idx) const
{
	assert(idx >= 0 && idx < (int)m_lensesAndRedshifts.size());
	return *(m_lensesAndRedshifts[idx].first.get());
}

inline double MultiPlaneContainerParams::getRedshift(int idx) const
{
	assert(idx >= 0 && idx < (int)m_lensesAndRedshifts.size());
	return m_lensesAndRedshifts[idx].second;
}

// Just a container to store multiple lenses at different redshifts
class GRALE_IMPORTEXPORT MultiPlaneContainer : public GravitationalLens
{
public:
	MultiPlaneContainer();
	~MultiPlaneContainer();

	bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const override;
	double getSurfaceMassDensity(Vector2D<double> theta) const override;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const override;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const override;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams) override;
private:
	static std::string s_onlyContainerError;
};

} // end namespace
