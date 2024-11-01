#include "mcmcparameters.h"
#include <array>

using namespace errut;
using namespace serut;
using namespace std;

namespace grale
{

MCMCParameters::MCMCParameters(double a, const std::string &samplesFileName,
	               size_t sampleGenerations, size_t burnInGenerations,
	               size_t annealGenerationsScale, double alpha0, double alphaMax)
	: EAParameters(EAParameters::MCMC),
	  m_a(a),
	  m_fileName(samplesFileName),
	  m_sampleGenerations(sampleGenerations),
	  m_burnInGenerations(burnInGenerations),
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
	array<int32_t,3> gens;
	if (!si.readInt32s(gens.data(), gens.size()) ||
	    !si.readDoubles(reals.data(), reals.size()) ||
		!si.readString(m_fileName))
		return si.getErrorString();

	m_a = reals[0];
	m_alpha0 = reals[1];
	m_alphaMax = reals[2];
	m_sampleGenerations = (size_t)gens[0];
	m_annealGenerations = (size_t)gens[1];
	m_burnInGenerations = (size_t)gens[2];
	return true;
}

bool_t MCMCParameters::writeInternal(SerializationInterface &si) const
{
	array<double,3> reals = { m_a, m_alpha0, m_alphaMax };
	array<int32_t,3> gens = { (int32_t)m_sampleGenerations, (int32_t)m_annealGenerations, (int32_t)m_burnInGenerations };

	if (!si.writeInt32s(gens.data(), gens.size()) ||
	    !si.writeDoubles(reals.data(), reals.size()) ||
		!si.writeString(m_fileName))
		return si.getErrorString();
	return true;
}


} // end namespace
