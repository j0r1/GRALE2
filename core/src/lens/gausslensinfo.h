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
 * \file gausslensinfo.h
 */

#ifndef GRALE_GAUSSLENSINFO_H

#define GRALE_GAUSSLENSINFO_H

#include "graleconfig.h"
#include "vector2d.h"

namespace grale
{
	
/** Information about a specific Gauss distribution, used by the multiple
 *  Gauss lens and the lens inversion routine.
 */
class GRALE_IMPORTEXPORT GaussLensInfo
{
public:
	GaussLensInfo()									{ m_mass = 0; m_angularWidth = 0; }

	/** Describe a Gauss mass distribution. */
	/** Describe a Gauss mass distribution.
	 *  \param mass The total mass of the distribution.
	 *  \param angularwidth The width of the distribution.
	 *  \param angularpos Position of the center of the mass distribution.
	 */
	GaussLensInfo(double mass, double angularWidth, Vector2D<double> angularPosition) { m_mass = mass; m_angularWidth = angularWidth; m_angularPosition = angularPosition; }

	/** Returns the currently set mass. */
	double getMass() const								{ return m_mass; }

	/** Returns the currently set width. */
	double getAngularWidth() const							{ return m_angularWidth; }

	/** Returns the currently set position of the center of the mass distribution. */
	Vector2D<double> getAngularPosition() const					{ return m_angularPosition; }

	/** Sets the current mass. */
	void setMass(double m)								{ m_mass = m; }
private:
	double m_mass, m_angularWidth;
	Vector2D<double> m_angularPosition;
};

} // end namespace

#endif // GRALE_GAUSSLENSINFO_H

