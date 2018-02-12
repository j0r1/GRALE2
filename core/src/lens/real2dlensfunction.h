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
 * \file real2dlensfunction.h
 */

#ifndef GRALE_REAL2DLENSFUNCTION_H

#define GRALE_REAL2DLENSFUNCTION_H

#include "graleconfig.h"
#include "real2dfunction.h"
#include "gravitationallens.h"

namespace grale
{

/** This class implements a 2D function which can be used to evaluate the mass density of a 
 *  lens at a specific point.
 */
class Real2DLensFunction : public Real2DFunction
{
public:	
	/** Initializes the function.
	 *  Initializes the function. Care must me taken that the lens \c l is valid
	 *  as long as the function is used.
	 *  \param l The lens of which the mass density should be calculated.
	 */
	Real2DLensFunction(GravitationalLens &l) : m_lens(l) 				{ }

	~Real2DLensFunction() 								{ }
	double operator()(Vector2D<double> theta) const					{ return m_lens.getSurfaceMassDensity(theta); }
private:
	GravitationalLens &m_lens;
};

} // end namespace

#endif // GRALE_REAL2DLENSFUNCTION_H

