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

/**
 * \file real1dfunctionintegrator.h
 */

#ifndef GRALE_REAL1DFUNCTIONINTEGRATOR_H

#define GRALE_REAL1DFUNCTIONINTEGRATOR_H

#include "graleconfig.h"
#include "real1dfunction.h"
#include <gsl/gsl_integration.h>
#include <vector>

namespace grale
{

/** Base class for an object which can integrate a function over some range. */
class GRALE_IMPORTEXPORT Real1DFunctionIntegrator
{
public:
	Real1DFunctionIntegrator(double absError = 0, double relError = 1e-5, int limit = 512);
	~Real1DFunctionIntegrator();

	double integrate(Real1DFunction &f, double x1, double x2);
private:
	static double staticIntegrationFunction(double x, void *pParams)
	{
		Real1DFunctionIntegrator *pInstance = (Real1DFunctionIntegrator *)pParams;

		return (*(pInstance->m_pFunction))(x);
	}

	gsl_integration_workspace *m_pWorkSpace;
	Real1DFunction *m_pFunction;
	double m_absError, m_relError;
	int m_limit;
};

} // end namespace

#endif // GRALE_REAL1DFUNCTIONINTEGRATOR_H

