#include "lensgaconvergenceparameters.h"

using namespace std;
using namespace serut;

namespace grale
{

LensGAConvergenceParameters::LensGAConvergenceParameters()
	: m_maxGenerations(16384),
	  m_convHistSize(250),
	  m_convFactor(0.1)
{
}

LensGAConvergenceParameters::~LensGAConvergenceParameters()
{
}

bool LensGAConvergenceParameters::write(SerializationInterface &si) const
{
	int32_t sizes[2] = { (int32_t)m_maxGenerations, (int32_t)m_convHistSize };
	if (!si.writeInt32s(sizes, 2))
	{
		setErrorString("Unable to write integer values: " + si.getErrorString());
		return false;
	}
	if (!si.writeDouble(m_convFactor))
	{
		setErrorString("Unable to write convergence factor: " + si.getErrorString());
		return false;
	}
	return true;
}

bool LensGAConvergenceParameters::read(SerializationInterface &si)
{
	int32_t sizes[2];
	if (!si.readInt32s(sizes, 2))
	{
		setErrorString("Unable to read integer values: " + si.getErrorString());
		return false;
	}
	for (auto s : sizes)
	{
		if (s < 0)
		{
			setErrorString("A negative size was read");
			return false;
		}
	}
	m_maxGenerations = (size_t)sizes[0];
	m_convHistSize = (size_t)sizes[1];

	if (!si.readDouble(&m_convFactor))
	{
		setErrorString("Error reading convergence factor: " + si.getErrorString());
		return false;
	}
	return true;
}

}

