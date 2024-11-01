#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

class MCMCParameters : public EAParameters
{
public:
	// For burnin, just don't output the samples
	// Annealing will calculate prob^alpha where alpha = alpha0 * (1/alpha0)^(t/T)
	MCMCParameters(double a = 2.0, const std::string &samplesFileName = "",
	               size_t sampleGenerations = 0,
	               size_t annealGenerationsScale = 0, double alpha0 = 0.1, double alphaMax = 1.0);
	~MCMCParameters();

	double getGoodmanWeare_a() const { return m_a; }
	const std::string &getSamplesFilename() const { return m_fileName; }
	size_t getSampleGenerations() const { return m_sampleGenerations; }
	size_t getAnnealGenerationsTimeScale() const { return m_annealGenerations; } // 0 means no annealing
	double getAnnealAlpha0() const { return m_alpha0; }
	double getAnnealAlphaMax() const { return m_alphaMax; }
protected:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_a;
	std::string m_fileName;
	size_t m_sampleGenerations;
	size_t m_annealGenerations;
	double m_alpha0, m_alphaMax;
};

} // end namespace
