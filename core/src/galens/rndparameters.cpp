#include "rndparameters.h"

using namespace errut;
using namespace serut;

namespace grale
{

RNDParameters::RNDParameters(double scale) : EAParameters(EAParameters::RND), m_scale(scale)
{
}

RNDParameters::~RNDParameters()
{
}

bool_t RNDParameters::readInternal(SerializationInterface &si)
{
	if (!si.readDouble(&m_scale))
		return si.getErrorString();
	return true;
}

bool_t RNDParameters::writeInternal(SerializationInterface &si) const
{
	if (!si.writeDouble(m_scale))
		return si.getErrorString();
	return true;
}

}
