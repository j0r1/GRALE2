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
	void setConvergenceFactor(double convFactor)							{ m_convFactor = convFactor; }

	size_t getMaximumNumberOfGenerations() const							{ return m_maxGenerations; }
	size_t getConvergenceHistorySize() const								{ return m_convHistSize; }
	double getConvergenceFactor() const										{ return m_convFactor; }

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;
private:
	size_t m_maxGenerations;
	size_t m_convHistSize;
	double m_convFactor;
};

}
