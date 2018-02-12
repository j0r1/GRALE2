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
 * \file rectangle2d.h
 */

#ifndef GRALE_RECTANGLE2D_H

#define GRALE_RECTANGLE2D_H

#include "graleconfig.h"
#include <vector>

namespace grale
{

/** Class which can hold the corner points of a rectangle. */
template <class T>
class Rectangle2D
{
public:
	Rectangle2D()										{ }
	~Rectangle2D()										{ } 

	void init(const Vector2D<T> points[4])							{ m_corners[0] = points[0]; m_corners[1] = points[1]; m_corners[2] = points[2]; m_corners[3] = points[3]; }
	const Vector2D<T> *getPoints() const							{ return m_corners; }
private:
	Vector2D<T> m_corners[4];
};
	
}

#endif // GRALE_RECTANGLE2D_H

