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

JADEParameters::JADEParameters(
		double p, double c,
		bool useArchive,
		double initMuF,
		double initMuCR
		)
	: EAParameters(EAParameters::JADE),
	  m_p(p), m_c(c),
	  m_useArchive(useArchive),
	  m_initMuF(initMuF),
	  m_initMuCR(initMuCR)
{
}

JADEParameters::~JADEParameters()
{
}

bool_t JADEParameters::readInternal(SerializationInterface &si)
{
	int32_t archInt;
	double params[4];

	if (!si.readInt32(&archInt) || !si.readDoubles(params, 4))
		return si.getErrorString();

	m_useArchive = (archInt == 0)?false:true;
	m_p = params[0];
	m_c = params[1];
	m_initMuF = params[2];
	m_initMuCR = params[3];
	return true;
}

bool_t JADEParameters::writeInternal(SerializationInterface &si) const
{
	int32_t archInt = (m_useArchive)?1:0;
	double params[4] = { m_p, m_c, m_initMuF, m_initMuCR };
	if (!si.writeInt32(archInt) || !si.writeDoubles(params, 4))
		return si.getErrorString();
	return true;
}

} // end namespace
