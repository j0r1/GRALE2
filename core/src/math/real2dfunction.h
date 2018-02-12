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
 * \file real2dfunction.h
 */

#ifndef GRALE_REAL2DFUNCTION_H

#define GRALE_REAL2DFUNCTION_H

#include "graleconfig.h"
#include "vector2d.h"

namespace grale
{

/** Base class for a 2D function. */
class GRALE_IMPORTEXPORT Real2DFunction
{
public:
	Real2DFunction()								{ }
	virtual ~Real2DFunction()							{ }
	virtual double operator()(Vector2D<double> v) const = 0;
};

}

#endif // GRALE_REAL2DFUNCTION_H

