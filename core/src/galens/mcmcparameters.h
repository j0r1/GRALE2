#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

class GeneralMCMCParameters : public EAParameters
{
protected:
	// For burnin, just don't output the samples
	// Annealing will calculate prob^alpha where alpha = alpha0 * (1/alpha0)^(t/T)
	GeneralMCMCParameters(ParameterType t, const std::string &samplesFileName = "",
				   const std::string &logProbFn = "",
	               size_t sampleGenerations = 0,
				   size_t burninGenerations = 0, // this many will not be logged to the sample file
	               size_t annealGenerationsScale = 0, double alpha0 = 0.1, double alphaMax = 1.0);
public:
	~GeneralMCMCParameters();

	const std::string getSamplesFilename() const { return m_fileName; }
	const std::string getLogProbFilename() const { return m_logProbFn; }
	// Note: these generations will be multiplied by two when starting the
	// actual algorithm, since two rounds are needed to advance all goodman-weare
	// walkers
	size_t getSampleGenerations() const { return m_sampleGenerations; }
	size_t getBurnInGenerations() const { return m_burnInGenerations; }
	size_t getAnnealGenerationsTimeScale() const { return m_annealGenerations; } // 0 means no annealing
	double getAnnealAlpha0() const { return m_alpha0; }
	double getAnnealAlphaMax() const { return m_alphaMax; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	std::string m_fileName, m_logProbFn;
	size_t m_sampleGenerations;
	size_t m_burnInGenerations;
	size_t m_annealGenerations;
	double m_alpha0, m_alphaMax;
};

// TODO: rename this to be specific for goodman-weare
class MCMCParameters : public GeneralMCMCParameters
{
public:
	MCMCParameters(double a = 2.0, const std::string &samplesFileName = "",
			       const std::string &logProbFn = "",
	               size_t sampleGenerations = 0,
				   size_t burninGenerations = 0, // this many will not be logged to the sample file
	               size_t annealGenerationsScale = 0, double alpha0 = 0.1, double alphaMax = 1.0);
	~MCMCParameters();

	double getGoodmanWeare_a() const { return m_a; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_a;
};

class MetropolisHastingsMCMCParameters : public GeneralMCMCParameters
{
public:
	MetropolisHastingsMCMCParameters(const std::vector<double> &stepScales = {},
				   const std::string &samplesFileName = "",
				   const std::string &logProbFn = "",
	               size_t sampleGenerations = 0,
				   size_t burninGenerations = 0, // this many will not be logged to the sample file
	               size_t annealGenerationsScale = 0, double alpha0 = 0.1, double alphaMax = 1.0);
	~MetropolisHastingsMCMCParameters();

	const std::vector<double> &getStepScales() const { return m_stepScales; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	std::vector<double> m_stepScales;
};

} // end namespace
