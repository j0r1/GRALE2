#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

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
