#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

class DEParameters : public EAParameters
{
public:
	DEParameters(double F = 0.5, double CR = 0.5, bool needStrictlyBetter = true);
	~DEParameters();

	double getF() const { return m_F; }
	double getCR() const { return m_CR; }
	bool getNeedStrictlyBetter() const { return m_needStrictlyBetter; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_F, m_CR;
	bool m_needStrictlyBetter;
};

class JADEParameters : public EAParameters
{
public:
	JADEParameters(
		double p = 0.05, double c = 0.1,
		bool useArchive = true,
		double initMuF = 0.5,
		double initMuCR = 0.5,
		bool needStrictlyBetter = true
		);
	~JADEParameters();

	double getBestFraction_p() const { return m_p; }
	double getParameterUpdateFraction_c() const { return m_c; }
	bool useExternalArchive() const { return m_useArchive; }
	double getInitialMeanF() const { return m_initMuF; }
	double getInitialMeanCR() const { return m_initMuCR; }
	bool getNeedStrictlyBetter() const { return m_needStrictlyBetter; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_p, m_c;
	bool m_useArchive;
	double m_initMuF, m_initMuCR;
	bool m_needStrictlyBetter;
};

} // end namespace
