#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

class RNDParameters : public EAParameters
{
public:
	RNDParameters(double scale = 1e-5);
	~RNDParameters();

	double getScale() const { return m_scale; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_scale;
};

} // end namespace
