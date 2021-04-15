#include "lensgaconvergenceparameters.h"

using namespace std;
using namespace serut;

namespace grale
{

LensGAConvergenceParameters::LensGAConvergenceParameters()
	: m_maxGenerations(16384),
	  m_convHistSize(250)
{
}

LensGAConvergenceParameters::~LensGAConvergenceParameters()
{
}

bool LensGAConvergenceParameters::setConvergenceFactorsAndMutationSizes(const vector<double> &convFactors, const vector<double> &mutSizes)
{
	if (convFactors.size() != mutSizes.size())
	{
		setErrorString("Both arrays should have the same length");
		return false;
	}
	if (convFactors.size() < 1)
	{
		setErrorString("At least one entry needs to be provided");
		return false;
	}
	m_convFactors = convFactors;
	m_convMutSizes = mutSizes;
	return true;
}

bool LensGAConvergenceParameters::write(SerializationInterface &si) const
{
	int32_t sizes[4] = { (int32_t)m_maxGenerations, (int32_t)m_convHistSize,
		                 (int32_t)m_convFactors.size(), (int32_t)m_convMutSizes.size() };
	if (!si.writeInt32s(sizes, 4))
	{
		setErrorString("Unable to write integer values: " + si.getErrorString());
		return false;
	}
	if (!si.writeDoubles(m_convFactors) || !si.writeDoubles(m_convMutSizes))
	{
		setErrorString("Unable to write convergence or mutation sizes: " + si.getErrorString());
		return false;
	}
	return true;
}

bool LensGAConvergenceParameters::read(SerializationInterface &si)
{
	int32_t sizes[4];
	if (!si.readInt32s(sizes, 4))
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
	m_convFactors.resize((size_t)sizes[2]);
	m_convMutSizes.resize((size_t)sizes[3]);

	if (!si.readDoubles(m_convFactors) || !si.readDoubles(m_convMutSizes))
	{
		setErrorString("Error reading convergence factors or mutation sizes: " + si.getErrorString());
		return false;
	}
	return true;
}

}
