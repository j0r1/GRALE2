#include "lensgamultipopulationparameters.h"

namespace grale
{

LensGAMultiPopulationParameters::LensGAMultiPopulationParameters()
	: m_numPop(2),
	  m_gracePeriod(250),
	  m_fraction(0.05),
	  m_iterations(1)
{
}

LensGAMultiPopulationParameters::~LensGAMultiPopulationParameters()
{
}

bool LensGAMultiPopulationParameters::read(serut::SerializationInterface &si)
{
	int32_t x[3];
	double f;

	if (!si.readInt32s(x, 3) || !si.readDouble(&f))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_numPop = (size_t)x[0];
	m_gracePeriod = (size_t)x[1];
	m_iterations = (size_t)x[2];
	m_fraction = f;

	return true;
}

bool LensGAMultiPopulationParameters::write(serut::SerializationInterface &si) const
{

	int32_t x[3] = { (int32_t)m_numPop, (int32_t)m_gracePeriod, (int32_t)m_iterations };
	if (!si.writeInt32s(x, 3) || !si.writeDouble(m_fraction))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

}
