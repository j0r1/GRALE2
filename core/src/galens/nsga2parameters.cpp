#include "nsga2parameters.h"

using namespace errut;
using namespace serut;

namespace grale
{

NSGA2Parameters::NSGA2Parameters(double smallMutationSize)
	: EAParameters(NSGA2), m_smallMutSize(smallMutationSize)
{
}

NSGA2Parameters::~NSGA2Parameters()
{
}

bool_t NSGA2Parameters::readInternal(serut::SerializationInterface &si)
{
	if (!si.readDouble(&m_smallMutSize))
		return si.getErrorString();
	return true;
}

bool_t NSGA2Parameters::writeInternal(serut::SerializationInterface &si) const
{
	if (!si.writeDouble(m_smallMutSize))
		return si.getErrorString();
	return true;
}

}
