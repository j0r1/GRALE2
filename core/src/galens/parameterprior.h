#pragma once

#include "graleconfig.h"
#include <serut/serializationinterface.h>
#include <errut/booltype.h>
#include <memory>

namespace grale
{

class ParameterPrior
{
public:
	enum PriorType { Invalid = 0, None, Gaussian };

	virtual ~ParameterPrior() { }

	static errut::bool_t read(serut::SerializationInterface &si, std::unique_ptr<ParameterPrior> &prior);
	errut::bool_t write(serut::SerializationInterface &si) const;

	virtual float getNegativeLogProb(float x) const = 0;
protected:
	ParameterPrior(PriorType t) : m_type(t) { }

	virtual errut::bool_t readInternal(serut::SerializationInterface &si) = 0;
	virtual errut::bool_t writeInternal(serut::SerializationInterface &si) const = 0;
private:
	PriorType m_type = Invalid;
};

class NoParameterPrior : public ParameterPrior
{
public:
	NoParameterPrior() : ParameterPrior(None) { }
	~NoParameterPrior() { }
private:
	float getNegativeLogProb(float x) const override { return 0; }

	errut::bool_t readInternal(serut::SerializationInterface &si) override { return true; }
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override { return true; }
};

class GaussianParameterPrior : public ParameterPrior
{
public:
	GaussianParameterPrior(float mu = 0, float sigma = 1) : ParameterPrior(Gaussian), m_mu(mu), m_sigma(sigma) { }
	~GaussianParameterPrior() { }
private:
	float getNegativeLogProb(float x) const override;

	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;

	float m_mu = 0, m_sigma = 1;
};

} // end namespace
