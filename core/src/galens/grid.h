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

#ifndef GRALE_GRID_H

#define GRALE_GRID_H

#include "graleconfig.h"
#include "real2dderivablefunction.h"
#include <errut/errorbase.h>
#include <list>

namespace grale
{

class GridSquare
{
public:
	GridSquare(Vector2D<double> center, double size)					{ m_center = center; m_size = size; }
	~GridSquare()										{ }
	double getSize() const									{ return m_size; }
	Vector2D<double> getCenter() const							{ return m_center; }
	Vector2D<double> getBottomLeft() const							{ return Vector2D<double>(m_center.getX()-m_size/2.0,m_center.getY()-m_size/2.0); }
	Vector2D<double> getTopRight() const							{ return Vector2D<double>(m_center.getX()+m_size/2.0,m_center.getY()+m_size/2.0); }
private:
	Vector2D<double> m_center;
	double m_size;
};

} // end namespace

#endif // GRALE_GRID_H

