#include "graleconfig.h"
#include "parameterprior.h"
#include <vector>

using namespace std;
using namespace errut;

namespace grale
{

bool_t ParameterPrior::read(serut::SerializationInterface &si, unique_ptr<ParameterPrior> &prior)
{
	int32_t tpe;
	if (!si.readInt32(&tpe))
		return "Can't read prior type: " + si.getErrorString();

	unique_ptr<ParameterPrior> newPrior;
	if (tpe == (int)ParameterPrior::Gaussian)
		newPrior = make_unique<GaussianParameterPrior>();
	else
		return "Invalid prior type " + to_string(tpe);

	bool_t r = newPrior->readInternal(si);
	if (!r)
		return r;
	prior = move(newPrior);
	return true;
}

bool_t ParameterPrior::write(serut::SerializationInterface &si) const
{
	if (!si.writeInt32((int)m_type))
		return "Can't write prior type: " + si.getErrorString();

	return writeInternal(si);
}

// Gaussian

float GaussianParameterPrior::getNegativeLogProb(float x) const
{
	float diff = (x-m_mu)/m_sigma;
	return 0.5*diff*diff;
}

bool_t GaussianParameterPrior::readInternal(serut::SerializationInterface &si)
{
	vector<float> vals(2);
	if (!si.readFloats(vals))
		return si.getErrorString();
	m_mu = vals[0];
	m_sigma = vals[1];
	if (m_sigma <= 0)
		return "Read invalid sigma: " + to_string(m_sigma);
	return true;
}

bool_t GaussianParameterPrior::writeInternal(serut::SerializationInterface &si) const
{
	vector<float> vals = { m_mu, m_sigma };
	if (!si.writeFloats(vals))
		return si.getErrorString();
	return true;
}


} // end namespace


