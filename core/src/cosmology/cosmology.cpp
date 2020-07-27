#include "cosmology.h"
#include "constants.h"
#include <limits>
#include <cmath>
#include <gsl/gsl_integration.h>

using namespace std;

namespace grale
{

Cosmology::Cosmology(double h, double Wm, double Wr, double Wv, double w)
{
    m_h = h;
    m_Wm = Wm;
    m_Wr = Wr;
    m_Wv = Wv;
    m_w = w;
}

Cosmology::~Cosmology()
{
    
}

struct CosmParams
{
	CosmParams(double _h, double _Wm, double _Wr, double _Wv, double _w)
	{
		h = _h;
		Wm = _Wm;
		Wr = _Wr;
		Wv = _Wv;
		w = _w;

		Wk = 1.0-Wm-Wr-Wv;
		if (std::abs(Wk) < 1.0e-7)
		{
			Wk = 0;
			Wv = 1.0-Wm-Wr;
		}
	}

	double h, Wm, Wr, Wv, Wk, w;
};

double Cosmology::getAngularDiameterDistance(double z1, double z2) const
{
    if (z1 < 0 || z2 < 0)
    {
        setErrorString("Redshifts must be positive");
        return numeric_limits<double>::quiet_NaN();
    }

    if (z1 > z2)
        std::swap(z1, z2);

   	CosmParams cosm(m_h, m_Wm, m_Wr, m_Wv, m_w);
	gsl_function F;
	F.function = integrationFunction;
	F.params = &cosm;

	// Check better for errors?
	double result = 0, abserr = 0;
	size_t neval = 0;

    int status = gsl_integration_qng(&F, 1.0/(1.0+z2), 1.0/(1.0+z1), 0, 1e-7, &result, &abserr, &neval);
	if (status)
    {
        setErrorString("Error performing gsl_integration_qng, code " + to_string(status));
        return numeric_limits<double>::quiet_NaN();
    }

	double T_H = DIST_MPC/100000.0;
	double commonFactor = 1.0/(1.0+z2) * SPEED_C*T_H/m_h;

	if (cosm.Wk == 0)
		return commonFactor * result;
	
	if (cosm.Wk < 0)
		return commonFactor/std::sqrt(-cosm.Wk) * std::sin(std::sqrt(-cosm.Wk) * result);
	
	// m_Wk > 0
	return commonFactor/std::sqrt(cosm.Wk) * std::sinh(std::sqrt(cosm.Wk) * result);
}

double Cosmology::integrationFunction(double R, void *params)
{
	CosmParams *pInst = reinterpret_cast<CosmParams*>(params);
	double R2 = R*R;
	return 1.0/std::sqrt(pInst->Wm*R + pInst->Wr + pInst->Wv*std::pow(R, 1.0-3.0*pInst->w) + pInst->Wk*R2);
}

bool Cosmology::write(serut::SerializationInterface &si) const
{
	double values[5] = { m_h, m_Wm, m_Wr, m_Wv, m_w };
	if (!si.writeDoubles(values, 5))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	return true;
}

bool Cosmology::read(serut::SerializationInterface &si)
{
	double values[5];
	if (!si.readDoubles(values, 5))
	{
		setErrorString(si.getErrorString());
		return false;
	}
	m_h = values[0];
	m_Wm = values[1];
	m_Wr = values[2];
	m_Wv = values[3];
	m_w = values[4];
	return true;
}

} // end namespace

