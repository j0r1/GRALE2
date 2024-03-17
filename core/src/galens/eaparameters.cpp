#include "eaparameters.h"
#include "gaparameters.h"
#include "deparameters.h"

using namespace errut;
using namespace std;
using namespace serut;

namespace grale
{

EAParameters::EAParameters(ParameterType t) : m_type(t)
{
}

EAParameters::~EAParameters()
{
}

bool_t EAParameters::read(SerializationInterface &si, unique_ptr<EAParameters> &parameters)
{
	int32_t typeInt;
	if (!si.readInt32(&typeInt))
		return "Can't read EAParameters type: " + si.getErrorString();

	unique_ptr<EAParameters> params;
	switch((ParameterType)typeInt)
	{
	case GA:
		params = make_unique<GAParameters>();
		break;
	case JADE:
		params = make_unique<JADEParameters>();
		break;
	default:
		return "Can't interpret EAParameters type id " + to_string(typeInt);
	}

	bool_t r = params->readInternal(si);
	if (!r)
		return r;

	parameters = move(params);
	return true;
}

bool_t EAParameters::write(SerializationInterface &si) const
{
	int32_t typeInt = (int32_t)m_type;
	if (!si.writeInt32(typeInt))
		return "Can't write EAParameters type: " + si.getErrorString();

	return writeInternal(si);
}

} // end namespace
