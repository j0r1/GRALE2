#include "scalesearchparameters.h"
#include <sstream>

using namespace std;

namespace grale
{

ScaleSearchParameters::ScaleSearchParameters(float startFactor, float stopFactor, int numIt, int firstItSteps, int subseqItSteps)
{
	m_startFactor = startFactor;
	m_stopFactor = stopFactor;
	m_numIterations = numIt;
	m_firstIterationsSteps = firstItSteps;
	m_subseqIterationSteps = subseqItSteps;
}

ScaleSearchParameters::ScaleSearchParameters(bool wideSearch)
{
	if (!wideSearch)
	{
		m_numIterations = 5;
		m_firstIterationsSteps = 20;
		m_subseqIterationSteps = 20;
		m_startFactor = 0.2;
		m_stopFactor = 5.0;
	}
	else
	{
		m_numIterations = 11;
		m_firstIterationsSteps = 50;
		m_subseqIterationSteps = 5;
		m_startFactor = 0.01;
		m_stopFactor = 100.0;
	}
}

ScaleSearchParameters::ScaleSearchParameters() // No search
{
	m_numIterations = 0;
	m_firstIterationsSteps = 0;
	m_subseqIterationSteps = 0;
	m_startFactor = 1.0;
	m_stopFactor = 1.0;
}

ScaleSearchParameters::~ScaleSearchParameters()
{
}

bool ScaleSearchParameters::write(serut::SerializationInterface &si) const
{
	int32_t iVals[] = { m_numIterations, m_firstIterationsSteps, m_subseqIterationSteps };
	float fVals[] = { m_startFactor, m_stopFactor };

	if (!si.writeInt32s(iVals, 3) || !si.writeFloats(fVals, 2))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool ScaleSearchParameters::read(serut::SerializationInterface &si)
{
	int32_t iVals[3];
	float fVals[2];

	if (!si.readInt32s(iVals, 3) || !si.readFloats(fVals, 2))
	{
		setErrorString(si.getErrorString());
		return false;
	}

	m_numIterations = iVals[0];
	m_firstIterationsSteps = iVals[1];
	m_subseqIterationSteps = iVals[2];
	m_startFactor = fVals[0];
	m_stopFactor = fVals[1];
	return true;
}

std::string ScaleSearchParameters::toString() const
{
	stringstream ss;
	ss << "f0=" << m_startFactor << " f1=" << m_stopFactor
	   << " ni=" << m_numIterations << " steps1=" << m_firstIterationsSteps
	   << " steps2=" << m_subseqIterationSteps;
	return ss.str();
}

bool ScaleSearchParameters::operator==(const ScaleSearchParameters &src) const
{
	return (m_startFactor == src.m_startFactor) && (m_stopFactor == src.m_stopFactor) &&
	       (m_numIterations == src.m_numIterations) &&
		   (m_firstIterationsSteps == src.m_firstIterationsSteps) &&
		   (m_subseqIterationSteps == src.m_subseqIterationSteps);
}

} // end namespace