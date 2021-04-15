#pragma once

#include "graleconfig.h"
#include <errut/errorbase.h>
#include <serut/serializationinterface.h>

namespace grale
{

class LensGAConvergenceParameters : public errut::ErrorBase
{
public:
	LensGAConvergenceParameters();
	~LensGAConvergenceParameters();

	void setMaximumNumberOfGenerations(size_t s)							{ m_maxGenerations = s; }
	void setHistorySize(size_t s)											{ m_convHistSize = s; }
	bool setConvergenceFactorsAndMutationSizes(const std::vector<double> &convFactors, const std::vector<double> &mutSizes);

	size_t getMaximumNumberOfGenerations() const							{ return m_maxGenerations; }
	size_t getConvergenceHistorySize() const								{ return m_convHistSize; }

	const std::vector<double> getConvergenceFactors() const					{ return m_convFactors; }
	// Negative means large mutation
	const std::vector<double> getConvergenceSmallMutationSizes() const		{ return m_convMutSizes; }

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;
private:
	size_t m_maxGenerations;
	size_t m_convHistSize;
	std::vector<double> m_convMutSizes, m_convFactors;
};

}
