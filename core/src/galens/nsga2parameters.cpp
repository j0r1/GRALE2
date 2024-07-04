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

NSGA2DELikeCrossoverParameters::NSGA2DELikeCrossoverParameters(bool extraParent,
	float F, float CR)
	: EAParameters(NSGA2DELikeCrossover), m_extraParent(extraParent), m_F(F), m_CR(CR)
{
}

NSGA2DELikeCrossoverParameters::~NSGA2DELikeCrossoverParameters()
{
}

bool_t NSGA2DELikeCrossoverParameters::readInternal(serut::SerializationInterface &si)
{
	int32_t flag;
	if (!si.readInt32(&flag) || !si.readFloat(&m_F) || !si.readFloat(&m_CR))
		return si.getErrorString();
	m_extraParent = (flag == 0)?false:true;
	return true;
}

bool_t NSGA2DELikeCrossoverParameters::writeInternal(serut::SerializationInterface &si) const
{
	int32_t flag = (m_extraParent)?1:0;
	if (!si.writeInt32(flag) || !si.writeFloat(m_F) || !si.writeFloat(m_CR))
		return si.getErrorString();
	return true;
}

}
