#include "mcmcparameters.h"
#include <array>

using namespace errut;
using namespace serut;
using namespace std;

namespace grale
{

GeneralMCMCParameters::GeneralMCMCParameters(ParameterType t, const std::string &samplesFileName,
				   const std::string &logProbFn,
	               size_t sampleGenerations, size_t burnInGenerations,
	               size_t annealGenerationsScale, double alpha0, double alphaMax)
	: EAParameters(t),
	  m_fileName(samplesFileName),
	  m_logProbFn(logProbFn),
	  m_sampleGenerations(sampleGenerations),
	  m_burnInGenerations(burnInGenerations),
	  m_annealGenerations(annealGenerationsScale),
	  m_alpha0(alpha0),
	  m_alphaMax(alphaMax)
{
}

GeneralMCMCParameters::~GeneralMCMCParameters()
{
}

bool_t GeneralMCMCParameters::readInternal(SerializationInterface &si)
{
	array<double,2> reals;
	array<int32_t,3> gens;
	if (!si.readInt32s(gens.data(), gens.size()) ||
	    !si.readDoubles(reals.data(), reals.size()) ||
		!si.readString(m_fileName) || !si.readString(m_logProbFn) )
		return si.getErrorString();

	m_alpha0 = reals[0];
	m_alphaMax = reals[1];
	m_sampleGenerations = (size_t)gens[0];
	m_annealGenerations = (size_t)gens[1];
	m_burnInGenerations = (size_t)gens[2];
	return true;
}

bool_t GeneralMCMCParameters::writeInternal(SerializationInterface &si) const
{
	array<double,2> reals = { m_alpha0, m_alphaMax };
	array<int32_t,3> gens = { (int32_t)m_sampleGenerations, (int32_t)m_annealGenerations, (int32_t)m_burnInGenerations };

	if (!si.writeInt32s(gens.data(), gens.size()) ||
	    !si.writeDoubles(reals.data(), reals.size()) ||
		!si.writeString(m_fileName) || !si.writeString(m_logProbFn) )
		return si.getErrorString();
	return true;
}

MCMCParameters::MCMCParameters(double a, const std::string &samplesFileName,
				   const std::string &logProbFn,
	               size_t sampleGenerations, size_t burnInGenerations,
	               size_t annealGenerationsScale, double alpha0, double alphaMax)
	: GeneralMCMCParameters(EAParameters::MCMC, samplesFileName, logProbFn, sampleGenerations, burnInGenerations, annealGenerationsScale, alpha0, alphaMax),
	  m_a(a)
{
}

MCMCParameters::~MCMCParameters()
{
}

bool_t MCMCParameters::readInternal(SerializationInterface &si)
{
	bool_t r = GeneralMCMCParameters::readInternal(si);
	if (!r)
		return r;

	if (!si.readDouble(&m_a))
		return si.getErrorString();

	return true;
}

bool_t MCMCParameters::writeInternal(SerializationInterface &si) const
{
	bool_t r = GeneralMCMCParameters::writeInternal(si);
	if (!r)
		return r;

	if (!si.writeDouble(m_a))
		return si.getErrorString();

	return true;
}


MetropolisHastingsMCMCParameters::MetropolisHastingsMCMCParameters(const std::vector<double> &stepScales,
			   const std::string &samplesFileName, const std::string &logProbFn,
			   size_t sampleGenerations,
			   size_t burninGenerations,
			   size_t annealGenerationsScale, double alpha0, double alphaMax)
	: GeneralMCMCParameters(EAParameters::MetropolisHastingsMCMC, samplesFileName, logProbFn, sampleGenerations, burninGenerations, annealGenerationsScale, alpha0, alphaMax),
	  m_stepScales(stepScales)
{
}

MetropolisHastingsMCMCParameters::~MetropolisHastingsMCMCParameters()
{
}

errut::bool_t MetropolisHastingsMCMCParameters::readInternal(serut::SerializationInterface &si)
{
	bool_t r = GeneralMCMCParameters::readInternal(si);
	if (!r)
		return r;

	int32_t num;
	if (!si.readInt32(&num))
		return si.getErrorString();

	if (num < 0)
		return "Read invalid number of step scales: " + to_string(num);

	m_stepScales.resize(num);
	if (!si.readDoubles(m_stepScales))
		return si.getErrorString();

	return true;
}

errut::bool_t MetropolisHastingsMCMCParameters::writeInternal(serut::SerializationInterface &si) const
{
	bool_t r = GeneralMCMCParameters::writeInternal(si);
	if (!r)
		return r;

	if (!si.writeInt32((int32_t)m_stepScales.size()) ||
	    !si.writeDoubles(m_stepScales))
		return si.getErrorString();

	return true;
}

} // end namespace
