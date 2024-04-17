#include "deparameters.h"

using namespace errut;
using namespace serut;

namespace grale
{

DEParameters::DEParameters(double F, double CR, bool needStrictlyBetter)
	: EAParameters(EAParameters::DE), m_F(F), m_CR(CR), m_needStrictlyBetter(needStrictlyBetter)
{
}

DEParameters::~DEParameters()
{
}

bool_t DEParameters::readInternal(SerializationInterface &si)
{
	int32_t betterInt;
	if (!si.readDouble(&m_F) || !si.readDouble(&m_CR) || !si.readInt32(&betterInt))
		return si.getErrorString();

	m_needStrictlyBetter = (betterInt == 0)?false:true;

	return true;
}

bool_t DEParameters::writeInternal(SerializationInterface &si) const
{
	int32_t betterInt = (m_needStrictlyBetter)?1:0;

	if (!si.writeDouble(m_F) || !si.writeDouble(m_CR) || !si.writeInt32(betterInt))
		return si.getErrorString();
	return true;
}

JADEParameters::JADEParameters(
		double p, double c,
		bool useArchive,
		double initMuF,
		double initMuCR,
		bool needStrictlyBetter
		)
	: EAParameters(EAParameters::JADE),
	  m_p(p), m_c(c),
	  m_useArchive(useArchive),
	  m_initMuF(initMuF),
	  m_initMuCR(initMuCR),
	  m_needStrictlyBetter(needStrictlyBetter)
{
}

JADEParameters::~JADEParameters()
{
}

bool_t JADEParameters::readInternal(SerializationInterface &si)
{
	int32_t archInt, betterInt;
	double params[4];

	if (!si.readInt32(&archInt) || !si.readDoubles(params, 4) || !si.readInt32(&betterInt))
		return si.getErrorString();

	m_useArchive = (archInt == 0)?false:true;
	m_p = params[0];
	m_c = params[1];
	m_initMuF = params[2];
	m_initMuCR = params[3];
	m_needStrictlyBetter = (betterInt == 0)?false:true;
	return true;
}

bool_t JADEParameters::writeInternal(SerializationInterface &si) const
{
	int32_t archInt = (m_useArchive)?1:0;
	int32_t betterInt = (m_needStrictlyBetter)?1:0;
	double params[4] = { m_p, m_c, m_initMuF, m_initMuCR };
	if (!si.writeInt32(archInt) || !si.writeDoubles(params, 4) || !si.writeInt32(betterInt))
		return si.getErrorString();
	return true;
}

} // end namespace
