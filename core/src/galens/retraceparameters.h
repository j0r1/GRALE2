#pragma once

#include "graleconfig.h"
#include <serut/serializationinterface.h>
#include <errut/booltype.h>
#include <memory>

namespace grale
{

class TraceParameters
{
public:
	enum ParameterType { SingleStepNewton, MultiStepNewton, ExpandedMultiStepNewton };

	virtual ~TraceParameters();

	static errut::bool_t read(serut::SerializationInterface &si, std::unique_ptr<TraceParameters> &parameters);
	errut::bool_t write(serut::SerializationInterface &si) const;
protected:
	TraceParameters(ParameterType t);

	virtual errut::bool_t readInternal(serut::SerializationInterface &si) = 0;
	virtual errut::bool_t writeInternal(serut::SerializationInterface &si) const = 0;
private:
	ParameterType m_type;
};

class SingleStepNewtonTraceParams : public TraceParameters
{
public:
	SingleStepNewtonTraceParams() : TraceParameters(SingleStepNewton) { }
	~SingleStepNewtonTraceParams() { }
private:
	errut::bool_t readInternal(serut::SerializationInterface &si) override { return true; }
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override { return true; }
};

class MultiStepNewtonTraceParams : public TraceParameters
{
public:
	MultiStepNewtonTraceParams(size_t numEvaluations = 0) : TraceParameters(MultiStepNewton), m_numEvals(numEvaluations) { }
	~MultiStepNewtonTraceParams() { }

	size_t getNumberOfEvaluations() const { return m_numEvals; }
private:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;

	size_t m_numEvals;
};

class ExpandedMultiStepNewtonTraceParams : public TraceParameters
{
public:
	ExpandedMultiStepNewtonTraceParams(size_t numEvalsPerStartPosition = 0, size_t numMaxGridSteps = 0, double acceptThreshold = 0,
									   double gridSpacing = 0)
		: TraceParameters(ExpandedMultiStepNewton),
		  m_numEvalsPerStartPos(numEvalsPerStartPosition), m_numMaxGridSteps(numMaxGridSteps),
	      m_acceptThreshold(acceptThreshold), m_gridSpacing(gridSpacing) { }
	~ExpandedMultiStepNewtonTraceParams() { }

	size_t getNumberOfEvaluationsPerStartPosition() const { return m_numEvalsPerStartPos; }
	size_t getMaximumNumberOfGridSteps() const { return m_numMaxGridSteps; }
	double getAcceptanceThreshold() const { return m_acceptThreshold; }
	double getGridSpacing() const { return m_gridSpacing; }
private:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;

	size_t m_numEvalsPerStartPos, m_numMaxGridSteps;
	double m_acceptThreshold, m_gridSpacing;
};

} // end namespace
