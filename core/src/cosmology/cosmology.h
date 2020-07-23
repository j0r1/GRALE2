#pragma once

#include "graleconfig.h"
#include <errut/errorbase.h>

namespace grale
{

class Cosmology : public errut::ErrorBase
{
public:
    Cosmology(double h, double Wm, double Wr, double Wv, double w = -1.0);
    ~Cosmology();

    double getH() const                                                     { return m_h; }
    double getOmegaM() const                                                { return m_Wm; }
    double getOmegaR() const                                                { return m_Wr; }
    double getOmegaV() const                                                { return m_Wv; }
    double getW() const                                                     { return m_w; }

    double getAngularDiameterDistance(double z1, double z2 = 0);
private:
	static double integrationFunction(double R, void *params);

    double m_h, m_Wm, m_Wr, m_Wv, m_w;
};

} // end namespace