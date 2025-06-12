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
	enum ParameterType { NoTrace, SingleStepNewton, MultiStepNewton, ExpandedMultiStepNewton };

	virtual ~TraceParameters();

	static errut::bool_t read(serut::SerializationInterface &si, std::unique_ptr<TraceParameters> &parameters);
	errut::bool_t write(serut::SerializationInterface &si) const;

	virtual std::string getRetraceDescription() const = 0;
	virtual std::unique_ptr<TraceParameters> createScaledCopy(double angScale, double potScale) const = 0;
protected:
	TraceParameters(ParameterType t);

	virtual errut::bool_t readInternal(serut::SerializationInterface &si) = 0;
	virtual errut::bool_t writeInternal(serut::SerializationInterface &si) const = 0;
private:
	ParameterType m_type;
};

class NoTraceParameters : public TraceParameters
{
public:
	NoTraceParameters() : TraceParameters(NoTrace) { }
	~NoTraceParameters() { }

	std::string getRetraceDescription() const override { return "NoTrace"; }
	std::unique_ptr<TraceParameters> createScaledCopy(double angScale, double potScale) const override { return std::make_unique<NoTraceParameters>(); }
private:
	errut::bool_t readInternal(serut::SerializationInterface &si) override { return true; }
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override { return true; }
};

class SingleStepNewtonTraceParams : public TraceParameters
{
public:
	SingleStepNewtonTraceParams() : TraceParameters(SingleStepNewton) { }
	~SingleStepNewtonTraceParams() { }

	std::string getRetraceDescription() const override { return "SingleStepNewton"; }
	std::unique_ptr<TraceParameters> createScaledCopy(double angScale, double potScale) const override { return std::make_unique<SingleStepNewtonTraceParams>(); }
private:
	errut::bool_t readInternal(serut::SerializationInterface &si) override { return true; }
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override { return true; }
};

class MultiStepNewtonTraceParams : public TraceParameters
{
public:
	MultiStepNewtonTraceParams(size_t numEvaluations = 0) : TraceParameters(MultiStepNewton), m_numEvals(numEvaluations) { }
	~MultiStepNewtonTraceParams() { }

	std::string getRetraceDescription() const override { return "MultiStepNewton, numEvals = " + std::to_string(m_numEvals); }
	size_t getNumberOfEvaluations() const { return m_numEvals; }
	std::unique_ptr<TraceParameters> createScaledCopy(double angScale, double potScale) const override { return std::make_unique<MultiStepNewtonTraceParams>(m_numEvals); }
private:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;

	size_t m_numEvals;
};

class ExpandedMultiStepNewtonTraceParams : public TraceParameters
{
public:
	enum Layout { Invalid, FullGrid, Square, Diamond, EightNeighbours };

	ExpandedMultiStepNewtonTraceParams(Layout l = Invalid, size_t numEvalsPerStartPosition = 0, size_t numMaxGridSteps = 0, double acceptThreshold = 0,
									   double gridSpacing = 0)
		: TraceParameters(ExpandedMultiStepNewton), m_layout(l),
		  m_numEvalsPerStartPos(numEvalsPerStartPosition), m_numMaxGridSteps(numMaxGridSteps),
	      m_acceptThreshold(acceptThreshold), m_gridSpacing(gridSpacing), m_rescaled(false) { }
	~ExpandedMultiStepNewtonTraceParams() { }

	std::string getRetraceDescription() const override;
	Layout getLayout() const { return m_layout; }
	size_t getNumberOfEvaluationsPerStartPosition() const { return m_numEvalsPerStartPos; }
	size_t getMaximumNumberOfGridSteps() const { return m_numMaxGridSteps; }
	double getAcceptanceThreshold() const { return m_acceptThreshold; }
	double getGridSpacing() const { return m_gridSpacing; }

	std::unique_ptr<TraceParameters> createScaledCopy(double angScale, double potScale) const override
	{
		auto copy = std::make_unique<ExpandedMultiStepNewtonTraceParams>(m_layout, m_numEvalsPerStartPos, m_numMaxGridSteps, m_acceptThreshold/angScale, m_gridSpacing/angScale); 
		copy->m_rescaled = true;
		return copy;
	}
private:
	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;

	Layout m_layout;
	size_t m_numEvalsPerStartPos, m_numMaxGridSteps;
	double m_acceptThreshold, m_gridSpacing;
	bool m_rescaled;
};

} // end namespace
