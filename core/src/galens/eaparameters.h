#pragma once

#include "graleconfig.h"
#include <serut/serializationinterface.h>
#include <errut/booltype.h>
#include <memory>

namespace grale
{

// Base class for different parameter types
class EAParameters
{
public:
	enum ParameterType { TEST, GA, DE, JADE, RND };
	virtual ~EAParameters();

	static errut::bool_t read(serut::SerializationInterface &si, std::unique_ptr<EAParameters> &parameters);
	errut::bool_t write(serut::SerializationInterface &si) const;
protected:
	EAParameters(ParameterType t);

	virtual errut::bool_t readInternal(serut::SerializationInterface &si) = 0;
	virtual errut::bool_t writeInternal(serut::SerializationInterface &si) const = 0;
private:
	ParameterType m_type;
};

class EATestParameters : public EAParameters
{
public:
	EATestParameters() : EAParameters(EAParameters::TEST) { }
	~EATestParameters() { }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override { return true; }
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override { return true; }
};

} // end namespace
