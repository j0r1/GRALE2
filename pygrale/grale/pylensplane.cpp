#include "pylensplane.h"
#include "images.h"
#include <serut/vectorserializer.h>
#include <grale/gravitationallens.h>

using namespace grale;
using namespace std;
using namespace serut;

PyLensPlane::PyLensPlane(PyObject *pPyObject) : m_pPyObject(pPyObject)
{
}

PyLensPlane::~PyLensPlane()
{
}

void PyLensPlane::setFeedbackStatus(const std::string &msg)
{
	if (PyObject_HasAttrString(m_pPyObject, "_onFeedbackStatus"))
		cy_call_void_string_function(m_pPyObject, "_onFeedbackStatus", (char *)msg.c_str());
}

void PyLensPlane::setFeedbackPercentage(int pct)
{
	if (PyObject_HasAttrString(m_pPyObject, "_onFeedbackPercentage"))
		cy_call_void_int_function(m_pPyObject, "_onFeedbackPercentage", pct);
}

bool PyLensPlane::createDeflectionGridLens(const LensPlane *lp, std::vector<uint8_t> &data, string &errStr)
{
	if (!lp)
	{
		errStr = "No lensplane specified";
		return false;
	}

	unique_ptr<GravitationalLens> pLens = lp->createDeflectionGridLens();
	if (!pLens.get())
	{
		errStr = lp->getErrorString();
		return false;
	}

	VectorSerializer vs;

	if (!pLens->write(vs))
	{
		errStr = "Unable to serialize lens: " + pLens->getErrorString();
		return false;
	}

	data = vs.getBuffer();

	return true;
}

template <class Plane>
bool getLensBytesTemplate(const Plane *lp, std::vector<uint8_t> &data, std::string &errStr)
{
	if (!lp)
	{
		errStr = "No lensplane or image plane specified";
		return false;
	}

	const GravitationalLens *pLens = lp->getLens();
	if (!pLens)
	{
		errStr = lp->getErrorString();
		return false;
	}

	VectorSerializer vs;
	if (!pLens->write(vs))
	{
		errStr = "Unable to serialize lens: " + pLens->getErrorString();
		return false;
	}

	data = vs.getBuffer();
	return true;
}

bool PyLensPlane::getLensBytes(const LensPlane *lp, std::vector<uint8_t> &data, std::string &errStr)
{
	return getLensBytesTemplate<LensPlane>(lp, data, errStr);
}

bool PyLensPlane::getLensBytesIP(const ImagePlane *lp, std::vector<uint8_t> &data, std::string &errStr)
{
	return getLensBytesTemplate<ImagePlane>(lp, data, errStr);
}

