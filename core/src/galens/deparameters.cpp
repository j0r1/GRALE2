#include "deparameters.h"

using namespace errut;
using namespace serut;

namespace grale
{

DEParameters::DEParameters(double F, double CR)
	: EAParameters(EAParameters::DE), m_F(F), m_CR(CR)
{
}

DEParameters::~DEParameters()
{
}

bool_t DEParameters::readInternal(SerializationInterface &si)
{
	if (!si.readDouble(&m_F) || !si.readDouble(&m_CR))
		return si.getErrorString();

	return true;
}

bool_t DEParameters::writeInternal(SerializationInterface &si) const
{
	if (!si.writeDouble(m_F) || !si.writeDouble(m_CR))
		return si.getErrorString();
	return true;
}

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
