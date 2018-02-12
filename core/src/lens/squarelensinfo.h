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
 * \file squarelensinfo.h
 */

#ifndef GRALE_SQUARELENSINFO_H

#define GRALE_SQUARELENSINFO_H

#include "graleconfig.h"
#include "vector2d.h"

namespace grale
{
	
/** Information about a specific square distribution, used by the multiple
 *  square lens and the lens inversion routine.
 */
class GRALE_IMPORTEXPORT SquareLensInfo
{
public:
	SquareLensInfo()								{ SquareLensInfo::mass = 0; SquareLensInfo::angularwidth = 0; SquareLensInfo::angularpos = Vector2D<double>(0,0); }

	/** Describe a square shaped mass distribution. */
	/** Describe a square shaped mass distribution.
	 *  \param mass The total mass of the distribution.
	 *  \param angularwidth The width of the distribution.
	 *  \param angularpos Position of the center of the mass distribution.
	 */
	SquareLensInfo(double mass, double angularWidth, Vector2D<double> angularPosition) { SquareLensInfo::mass = mass; SquareLensInfo::angularwidth = angularWidth; SquareLensInfo::angularpos = angularPosition; }

	/** Returns the currently set mass. */
	double getMass() const								{ return mass; }

	/** Returns the currently set width. */
	double getAngularWidth() const							{ return angularwidth; }

	/** Returns the currently set position of the center of the mass distribution. */
	Vector2D<double> getAngularPosition() const					{ return angularpos; }

	/** Sets the current mass. */
	void setMass(double m)								{ mass = m; }
private:
	double mass,angularwidth;
	Vector2D<double> angularpos;
};

} // end namespace

#endif // GRALE_SQUARELENSINFO_H

