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
#include "real1dfunctionintegrator.h"

namespace grale
{

Real1DFunctionIntegrator::Real1DFunctionIntegrator(double absError, double relError, int limit)
{
	m_absError = absError;
	m_relError = relError;
	m_limit = limit;

	m_pWorkSpace = gsl_integration_workspace_alloc(m_limit);
}

Real1DFunctionIntegrator::~Real1DFunctionIntegrator()
{
	gsl_integration_workspace_free(m_pWorkSpace);
}

double Real1DFunctionIntegrator::integrate(Real1DFunction &f, double x1, double x2)
{
	m_pFunction = &f;

	gsl_function F;
	
	F.function = staticIntegrationFunction;
	F.params = this;

	double result = 0;
	double err = 0;

	gsl_integration_qags(&F, x1, x2, m_absError, m_relError, m_limit, m_pWorkSpace, &result, &err); 

	return result;
}

} // end namespace


