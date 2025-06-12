#include "retraceparameters.h"
#include "constants.h"

using namespace std;
using namespace errut;

namespace grale
{

TraceParameters::TraceParameters(ParameterType t) : m_type(t)
{
}

TraceParameters::~TraceParameters()
{
}

errut::bool_t TraceParameters::read(serut::SerializationInterface &si, std::unique_ptr<TraceParameters> &parameters)
{
	int32_t typeInt = 0;

	if (!si.readInt32(&typeInt))
		return "Unable to read trace parameter type: " + si.getErrorString();

	unique_ptr<TraceParameters> params;
	switch((ParameterType)typeInt)
	{
	case NoTrace:
		params = make_unique<NoTraceParameters>();
		break;
	case SingleStepNewton:
		params = make_unique<SingleStepNewtonTraceParams>();
		break;
	case MultiStepNewton:
		params = make_unique<MultiStepNewtonTraceParams>();
		break;
	case ExpandedMultiStepNewton:
		params = make_unique<ExpandedMultiStepNewtonTraceParams>();
		break;
	default:
		return "Unknown trace parameter type " + to_string(typeInt);
	}

	bool_t r = params->readInternal(si);
	if (!r)
		return r;

	parameters = std::move(params);
	return true;
}

errut::bool_t TraceParameters::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32((int32_t)m_type))
		return "Unable to write trace parameter type: " + si.getErrorString();

	return writeInternal(si);
}

errut::bool_t MultiStepNewtonTraceParams::readInternal(serut::SerializationInterface &si)
{
	int32_t num;
	if (!si.readInt32(&num))
		return "Can't read number of evaluations: " + si.getErrorString();
	if (num <= 0)
		return "Invalid number of evaluations: " + to_string(num);

	m_numEvals = (size_t)num;
	return true;
}

errut::bool_t MultiStepNewtonTraceParams::writeInternal(serut::SerializationInterface &si) const
{
	int32_t num = (int32_t)m_numEvals;
	if (!si.writeInt32(num))
		return "Can't write number of evaluations: " + si.getErrorString();

	return true;
}

errut::bool_t ExpandedMultiStepNewtonTraceParams::readInternal(serut::SerializationInterface &si)
{
	vector<int32_t> intParams(3);
	vector<double> dblParams(2);

	if (!si.readInt32s(intParams) || !si.readDoubles(dblParams))
		return "Error reading integer or double parameters: " + si.getErrorString();

	for (int32_t x : intParams)
		if (x <= 0)
			return "Each integer parameter should be at least one, but detected " + to_string(x);
	
	m_layout = (Layout)intParams[0];
	m_numEvalsPerStartPos = (size_t)intParams[1];
	m_numMaxGridSteps = (size_t)intParams[2];

	for (double x : dblParams)
		if (x <= 0)
			return "Each floating point parameter should be strictly positive, but detected " + to_string(x);

	m_acceptThreshold = dblParams[0];
	m_gridSpacing = dblParams[1];
	return true;
}

errut::bool_t ExpandedMultiStepNewtonTraceParams::writeInternal(serut::SerializationInterface &si) const
{
	vector<int32_t> intParams { (int32_t)m_layout, (int32_t)m_numEvalsPerStartPos, (int32_t)m_numMaxGridSteps };
	vector<double> dblParams { m_acceptThreshold, m_gridSpacing };

	if (!si.writeInt32s(intParams) || !si.writeDoubles(dblParams))
		return "Error writing integer or double parameters: " + si.getErrorString();
	return true;
}

inline string getLayoutName(ExpandedMultiStepNewtonTraceParams::Layout l)
{
	if (l == ExpandedMultiStepNewtonTraceParams::FullGrid)
		return "FullGrid";
	if (l == ExpandedMultiStepNewtonTraceParams::Square)
		return "Square";
	if (l == ExpandedMultiStepNewtonTraceParams::Diamond)
		return "Diamond";
	if (l == ExpandedMultiStepNewtonTraceParams::EightNeighbours)
		return "EightNeighbours";
	return "Unknown";
}

std::string ExpandedMultiStepNewtonTraceParams::getRetraceDescription() const
{
	if (m_rescaled)
		return "ExpandedMultiStepNewton, numEvalsPerStartPos = " + to_string(m_numEvalsPerStartPos)
		   + ", numMaxGridSteps = " + to_string(m_numMaxGridSteps)
		   + ", acceptThreshold = " + to_string(m_acceptThreshold)
		   + ", gridSpacing = " + to_string(m_gridSpacing)
		   + ", layout = " + getLayoutName(m_layout);

	return "ExpandedMultiStepNewton, numEvalsPerStartPos = " + to_string(m_numEvalsPerStartPos)
		   + ", numMaxGridSteps = " + to_string(m_numMaxGridSteps)
		   + ", acceptThreshold = " + to_string(m_acceptThreshold/ANGLE_ARCSEC)
		   + " arcsec, gridSpacing = " + to_string(m_gridSpacing/ANGLE_ARCSEC) + " arcsec"
		   + ", layout = " + getLayoutName(m_layout);
}

} // end namespace
