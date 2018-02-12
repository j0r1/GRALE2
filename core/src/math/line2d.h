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
 * \file line2d.h
 */

#ifndef GRALE_LINE2D_H

#define GRALE_LINE2D_H

#include "graleconfig.h"
#include "vector2d.h"

namespace grale
{

/** Class which can represent a line. */
template<class T>
class Line2D
{
public:
	Line2D()										{ }
	Line2D(Vector2D<T> point, Vector2D<T> direction)					{ m_pt = point; m_dir = direction; }
	~Line2D()										{ }
	Vector2D<T> getIntersection(const Line2D<T> &l) const;
	Vector2D<T> getIntersection(const Line2D<T> &l, T &factor) const;
	T getIntersectionFactor(const Line2D<T> &l) const;
	Vector2D<T> getPoint(T factor) const							{ return m_pt + m_dir*factor; }
	Vector2D<T> getDirection() const							{ return m_dir; }
private:
	Vector2D<T> m_pt, m_dir;
};
	
template<class T>
inline Vector2D<T> Line2D<T>::getIntersection(const Line2D<T> &l) const
{
	Vector2D<T> diff = l.m_pt - m_pt;
	T coeff = (l.m_dir.getY()*diff.getX() - l.m_dir.getX()*diff.getY())/(m_dir.getX()*l.m_dir.getY() - m_dir.getY()*l.m_dir.getX());
	Vector2D<T> result = m_pt + m_dir*coeff;
	return result;
}

template<class T>
inline Vector2D<T> Line2D<T>::getIntersection(const Line2D<T> &l, T &factor) const
{
	Vector2D<T> diff = l.m_pt - m_pt;
	T coeff = (l.m_dir.getY()*diff.getX() - l.m_dir.getX()*diff.getY())/(m_dir.getX()*l.m_dir.getY() - m_dir.getY()*l.m_dir.getX());
	Vector2D<T> result = m_pt + m_dir*coeff;
	factor = coeff;
	return result;
}

template<class T>
inline T Line2D<T>::getIntersectionFactor(const Line2D<T> &l) const
{
	Vector2D<T> diff = l.m_pt-m_pt;
	T coeff = (l.m_dir.getY()*diff.getX() - l.m_dir.getX()*diff.getY())/(m_dir.getX()*l.m_dir.getY() - m_dir.getY()*l.m_dir.getX());
	return coeff;
}

} // end namespace

#endif // GRALE_LINE2D_H

