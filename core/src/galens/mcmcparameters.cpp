#include "mcmcparameters.h"
#include <array>

using namespace errut;
using namespace serut;
using namespace std;

namespace grale
{

MCMCParameters::MCMCParameters(double a, const std::string &samplesFileName,
	               size_t sampleGenerations,
	               size_t annealGenerationsScale, double alpha0, double alphaMax)
	: EAParameters(EAParameters::MCMC),
	  m_a(a),
	  m_fileName(samplesFileName),
	  m_sampleGenerations(sampleGenerations),
	  m_annealGenerations(annealGenerationsScale),
	  m_alpha0(alpha0),
	  m_alphaMax(alphaMax)
{
}

MCMCParameters::~MCMCParameters()
{
}

bool_t MCMCParameters::readInternal(SerializationInterface &si)
{
	array<double,3> reals;
	array<int32_t,2> gens;
	if (!si.readInt32s(gens.data(), gens.size()) ||
	    !si.readDoubles(reals.data(), reals.size()) ||
		!si.readString(m_fileName))
		return si.getErrorString();
	return true;
}

bool_t MCMCParameters::writeInternal(SerializationInterface &si) const
{
	array<double,3> reals = { m_a, m_alpha0, m_alphaMax };
	array<int32_t,2> gens = { (int32_t)m_sampleGenerations, (int32_t)m_annealGenerations };

	if (!si.writeInt32s(gens.data(), gens.size()) ||
	    !si.writeDoubles(reals.data(), reals.size()) ||
		!si.writeString(m_fileName))
		return si.getErrorString();
	return true;
}


} // end namespace
