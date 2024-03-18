#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

class DEParameters : public EAParameters
{
public:
	DEParameters(double F = 0.5, double CR = 0.5);
	~DEParameters();

	double getF() const { return m_F; }
	double getCR() const { return m_CR; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_F, m_CR;
};

class JADEParameters : public EAParameters
{
public:
	JADEParameters();
	~JADEParameters();
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
};

} // end namespace
