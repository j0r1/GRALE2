#include "deparameters.h"

using namespace errut;
using namespace serut;

namespace grale
{

JADEParameters::JADEParameters() : EAParameters(EAParameters::JADE)
{
}

JADEParameters::~JADEParameters()
{
}

bool_t JADEParameters::readInternal(SerializationInterface &si)
{
	return true;
}

bool_t JADEParameters::writeInternal(SerializationInterface &si) const
{
	return true;
}

} // end namespace
