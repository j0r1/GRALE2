#ifndef PYLENSPLANE_H

#define PYLENSPLANE_H

#include <Python.h>
#include <grale/lensplane.h>
#include <vector>

class PyLensPlane : public grale::LensPlane
{
public:
	PyLensPlane(PyObject *pPyObject);
	~PyLensPlane();

	static bool createDeflectionGridLens(const LensPlane *lp, std::vector<uint8_t> &data, std::string &errStr);
	static bool getLensBytes(const LensPlane *lp, std::vector<uint8_t> &data, std::string &errStr);
protected:
	void setFeedbackStatus(const std::string &msg);
	void setFeedbackPercentage(int pct);
private:
	PyObject *m_pPyObject;
};

#endif // PYLENSPLANE_H
