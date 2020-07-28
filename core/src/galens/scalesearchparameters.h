#pragma once

#include "graleconfig.h"
#include <errut/errorbase.h>
#include <serut/serializationinterface.h>

namespace grale
{

class GRALE_IMPORTEXPORT ScaleSearchParameters : public errut::ErrorBase
{
public:
	ScaleSearchParameters(float startFactor, float stopFactor, int numIt, int firstItSteps, int subseqItSteps);
	ScaleSearchParameters(bool wideSearch);
	ScaleSearchParameters(); // No search
	~ScaleSearchParameters();

	bool operator==(const ScaleSearchParameters &src) const;

	float getStartFactor() const 									{ return m_startFactor; }
	float getStopFactor() const 									{ return m_stopFactor; }
	int getNumberOfIterations() const 								{ return m_numIterations; }
	int getStepsOnFirstIteration() const 							{ return m_firstIterationsSteps; }
	int getStepsOnSubsequentIterations() const 						{ return m_subseqIterationSteps; }

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	std::string toString() const;
private:
	float m_startFactor, m_stopFactor;
	int m_firstIterationsSteps, m_subseqIterationSteps;
	int m_numIterations;
};

} // end namespace