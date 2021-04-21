/*

  This file is a part of GRALE, a library to facilitate the simulation
  and inversion of gravitational lenses.

  Copyright (C) 2008-2012 Jori Liesenborgs

  Contact: jori.liesenborgs@gmail.com
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
  
*/

#include "graleconfig.h"
#include "ellipticlens.h"
#include "real2dfunction.h"
#include "real1dfunction.h"
#include "real1dfunctionintegrator.h"
#include "circularlensprofile.h"
#include "constants.h"

#include <iostream>

namespace grale
{

class Cached2DFunction : public Real2DFunction
{
public:
	Cached2DFunction(Real2DFunction &f) : m_function(f)
	{
		m_hasCache = false;
	}

	double operator()(Vector2D<double> v) const
	{
		if (m_hasCache && m_cachePoint.getX() == v.getX() && m_cachePoint.getY() == v.getY())
		{
			return m_cacheValue;
		}

		double newValue = m_function(v);

		m_cacheValue = newValue;
		m_cachePoint = v;
		m_hasCache = true;

		return newValue;
	}
private:
	Real2DFunction &m_function;
	mutable Vector2D<double> m_cachePoint;
	mutable double m_cacheValue;
	mutable bool m_hasCache;
};

class EllipticIntegrand : public Real1DFunction
{
public:
	EllipticIntegrand(double q)
	{
		m_q = q;
		m_q2 = q*q;
		m_pCircularLensProfile = 0;
	}
	~EllipticIntegrand() 
	{ 
	}

	void setXY(double x, double y)
	{
		m_x = x;
		m_y = y;
		m_x2 = x*x;
		m_y2 = y*y;
	}

	void getXiAndDenom(double u, double &xi, double &denom) const
	{
		double d = 1.0-(1.0-m_q2)*u;
		double xi2 = u*(m_x2+m_y2/d);
		denom = d;
		xi = SQRT(xi2);
	}

	CircularLensProfile &getProfile() const
	{
		return *m_pCircularLensProfile;
	}

	void setProfile(CircularLensProfile *pProfile)
	{
		m_pCircularLensProfile = pProfile;
	}
private:
	double m_x, m_y, m_x2, m_y2, m_q, m_q2;
	CircularLensProfile *m_pCircularLensProfile;
};

class IntegrandI : public EllipticIntegrand
{
public:
	IntegrandI(double q) : EllipticIntegrand(q) { }
	~IntegrandI() { }

	double operator()(double u) const
	{
		double xi, denom;
		getXiAndDenom(u, xi, denom);

		double M = getProfile().getMassInside(xi);
		
		return (M/u)*(1.0/SQRT(denom));
	}
};

class IntegrandJ0 : public EllipticIntegrand
{
public:
	IntegrandJ0(double q) : EllipticIntegrand(q) { }
	~IntegrandJ0() { }

	double operator()(double u) const
	{
		double xi, denom;
		getXiAndDenom(u, xi, denom);

		double Sigma = getProfile().getSurfaceMassDensity(xi);
		
		return Sigma/SQRT(denom);
	}
};

class IntegrandJ1 : public EllipticIntegrand
{
public:
	IntegrandJ1(double q) : EllipticIntegrand(q) { }
	~IntegrandJ1() { }

	double operator()(double u) const
	{
		double xi, denom;
		getXiAndDenom(u, xi, denom);

		double Sigma = getProfile().getSurfaceMassDensity(xi);
		
		return Sigma/(SQRT(denom)*denom);
	}
};

class IntegrandK0 : public EllipticIntegrand
{
public:
	IntegrandK0(double q) : EllipticIntegrand(q) { }
	~IntegrandK0() { }

	double operator()(double u) const
	{
		double xi, denom;
		getXiAndDenom(u, xi, denom);

		double SigmaDeriv = getProfile().getSurfaceMassDensityDerivativeOverTheta(xi);
		
		return SigmaDeriv*u/SQRT(denom);
	}
};

class IntegrandK1 : public EllipticIntegrand
{
public:
	IntegrandK1(double q) : EllipticIntegrand(q) { }
	~IntegrandK1() { }

	double operator()(double u) const
	{
		double xi, denom;
		getXiAndDenom(u, xi, denom);

		double SigmaDeriv = getProfile().getSurfaceMassDensityDerivativeOverTheta(xi);
		
		return SigmaDeriv*u/(SQRT(denom)*denom);
	}
};

class IntegrandK2 : public EllipticIntegrand
{
public:
	IntegrandK2(double q) : EllipticIntegrand(q) { }
	~IntegrandK2() { }

	double operator()(double u) const
	{
		double xi, denom;
		getXiAndDenom(u, xi, denom);

		double SigmaDeriv = getProfile().getSurfaceMassDensityDerivativeOverTheta(xi);
		
		return SigmaDeriv*u/(SQRT(denom)*denom*denom);
	}
};

class EllipticFunction : public Real2DFunction
{
public:
	EllipticFunction(EllipticIntegrand &integrand, double absError, double relError, int limit) : m_integrand(integrand), m_integrator(absError, relError, limit)
	{
	}

	~EllipticFunction()
	{
	}

	double operator()(Vector2D<double> v) const
	{
		m_integrand.setXY(v.getX(), v.getY());
		return m_integrator.integrate(m_integrand, 0, 1);
	}
private:
	EllipticIntegrand &m_integrand;
	mutable Real1DFunctionIntegrator m_integrator;
};

template<class T>
class EllFunction : public Real2DFunction
{
public:
	EllFunction(double q, double absError, double relError, int limit, CircularLensProfile *pProfile) : m_integrand(q), m_function(m_integrand, absError, relError, limit), m_cache(m_function)
	{
		m_integrand.setProfile(pProfile);
	}

	double operator()(Vector2D<double> v) const
	{
		return m_cache(v);
	}
private:
	T m_integrand;
	EllipticFunction m_function;
	Cached2DFunction m_cache;
};

EllipticLens::EllipticLens(GravitationalLens::LensType t) : GravitationalLens(t)
{
}

EllipticLens::~EllipticLens()
{
}

bool EllipticLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const
{
	//std::cerr << "EllipticLens::getAlphaVector" << std::endl;
	double J0 = (*m_pJ0)(theta);
	double J1 = (*m_pJ1)(theta);

	double factor = ((4.0*CONST_G*CONST_PI)/(SPEED_C*SPEED_C))*getLensDistance()*m_q;

	Vector2D<double> a(factor*theta.getX()*J0, factor*theta.getY()*J1);
	*pAlpha = a;
	return true;
}

bool EllipticLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	//std::cerr << "EllipticLens::getAlphaVectorDerivatives" << std::endl;
	double J0 = (*m_pJ0)(theta);
	double J1 = (*m_pJ1)(theta);
	double K0 = (*m_pK0)(theta);
	double K1 = (*m_pK1)(theta);
	double K2 = (*m_pK2)(theta);

	double factor = ((4.0*CONST_G*CONST_PI)/(SPEED_C*SPEED_C))*getLensDistance()*m_q;
	double x = theta.getX();
	double y = theta.getY();
	
	axx = factor*x*x*K0 + factor*J0;
	ayy = factor*y*y*K2 + factor*J1;
	axy = factor*x*y*K1;

	return true;
}

double EllipticLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double x = theta.getX();
	double y = theta.getY()/m_q;
	double r = SQRT(x*x+y*y);
	return m_pProfile->getSurfaceMassDensity(r);
}

bool EllipticLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const
{
	//std::cerr << "EllipticLens::getProjectedPotential" << std::endl;
	double I = (*m_pI)(theta);
	double phi = (D_ds/D_s)*4.0*CONST_G/(SPEED_C*SPEED_C*getLensDistance())*m_q*0.5*I;
	*pPotentialValue = phi;
	return true;
}

void EllipticLens::subInit(double q, CircularLensProfile *pProfile,
		           double absError,
			   double relError,
			   int limit)
{
	m_q = q;
	m_pProfile = pProfile;

	m_pI = std::make_unique<EllFunction<IntegrandI>>(q, absError, relError, limit, pProfile);
	m_pJ0 = std::make_unique<EllFunction<IntegrandJ0>>(q, absError, relError, limit, pProfile);
	m_pJ1 = std::make_unique<EllFunction<IntegrandJ1>>(q, absError, relError, limit, pProfile);
	m_pK0 = std::make_unique<EllFunction<IntegrandK0>>(q, absError, relError, limit, pProfile);
	m_pK1 = std::make_unique<EllFunction<IntegrandK1>>(q, absError, relError, limit, pProfile);
	m_pK2 = std::make_unique<EllFunction<IntegrandK2>>(q, absError, relError, limit, pProfile);
}

} // end namespace
