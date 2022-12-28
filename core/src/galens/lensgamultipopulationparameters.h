#pragma once

#include "graleconfig.h"
#include <errut/errorbase.h>
#include <serut/serializationinterface.h>

namespace grale
{

class LensGAMultiPopulationParameters : public errut::ErrorBase
{
public:
	LensGAMultiPopulationParameters();
	~LensGAMultiPopulationParameters();

	void setNumberOfPopulations(size_t s) { m_numPop = s; }
	size_t getNumberOfPopulations() const { return m_numPop; }

	void setNumberOfInitialGenerationsToSkip(size_t n) { m_gracePeriod = n; }
	size_t getNumberOfInitialGenerationsToSkip() const { return m_gracePeriod; }

	void setMigrationGenerationFraction(double x) { m_fraction = x; }
	double getMigrationGenerationFraction() const { return m_fraction; }

	void setNumberOfIndividualsToLeavePopulation(size_t n) { m_iterations = n; }
	size_t getNumberOfIndividualsToLeavePopulation() const { return m_iterations; }

	bool read(serut::SerializationInterface &si);
	bool write(serut::SerializationInterface &si) const;
private:
	size_t m_numPop;
	size_t m_gracePeriod;
	double m_fraction;
	size_t m_iterations;
};

}
