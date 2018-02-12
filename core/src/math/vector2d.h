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
 * \file vector2d.h
 */

#ifndef GRALE_VECTOR2D_H

#define GRALE_VECTOR2D_H

#include "graleconfig.h"
#include "mathfunctions.h"
#include <assert.h>

namespace grale
{

/** Represents a 2D vector. */
template <class T>
class Vector2D
{
public:
	Vector2D()											{ m_pos[0] = 0; m_pos[1] = 0; }
	Vector2D(T x, T y)									{ m_pos[0] = x; m_pos[1] = y; }
	T getX() const										{ return m_pos[0]; }
	T getY() const										{ return m_pos[1]; }
	T getLength() const							
	{ 
		T absX = std::abs(m_pos[0]);
		T absY = std::abs(m_pos[1]);

		if (absX > absY)
		{
			T tmp = absY/absX;
			return absX*std::sqrt((T)1.0+tmp*tmp);
		}
		// absX <= absY
		if (absY == 0) // => absx == 0
			return 0;
		T tmp = absX/absY;
		return absY*std::sqrt(tmp*tmp+(T)1.0); 
	}
	T getLengthSquared() const							{ return m_pos[0]*m_pos[0]+m_pos[1]*m_pos[1]; }
	T getComponent(int i) const							{ if (i == 1) return m_pos[0]; assert(i == 2); return m_pos[1]; }
	Vector2D &operator+=(const Vector2D v)				{ m_pos[0] += v.m_pos[0]; m_pos[1] += v.m_pos[1]; return *this; }
	Vector2D &operator-=(const Vector2D v)				{ m_pos[0] -= v.m_pos[0]; m_pos[1] -= v.m_pos[1]; return *this; }
	Vector2D &operator*=(T a)							{ m_pos[0] *= a; m_pos[1] *= a; return *this; }
	Vector2D &operator/=(T a)							{ m_pos[0] /= a; m_pos[1] /= a; return *this; }
	const T *getComponents() const						{ return m_pos; }
	T *getComponents()									{ return m_pos; }
private:
	T m_pos[2];
};

template <class T>
inline Vector2D<T> operator+(Vector2D<T> v1, Vector2D<T> v2)
{
	return Vector2D<T>(v1.getX()+v2.getX(), v1.getY()+v2.getY());
}

template <class T>
inline Vector2D<T> operator-(Vector2D<T> v1, Vector2D<T> v2)
{
	return Vector2D<T>(v1.getX()-v2.getX(), v1.getY()-v2.getY());
}

template <class T>
inline Vector2D<T> operator*(Vector2D<T> v, T m)
{
	return Vector2D<T>(v.getX()*m, v.getY()*m);
}

template <class T>
inline Vector2D<T> operator*(T m, Vector2D<T> v)
{
	return Vector2D<T>(v.getX()*m, v.getY()*m);
}

template <class T>
inline Vector2D<T> operator/(Vector2D<T> v, T d)
{
	return Vector2D<T>(v.getX()/d, v.getY()/d);
}

template <class T>
inline T operator*(Vector2D<T> v1, Vector2D<T> v2)
{
	return v1.getX()*v2.getX() + v1.getY()*v2.getY();
}

typedef Vector2D<double> Vector2Dd;
typedef Vector2D<float> Vector2Df;

/** Represents a 2D vector with integer endpoint coordinates. */
class IntVector2D
{
public:
	IntVector2D()									{ m_x = 0; m_y = 0; }
	IntVector2D(int x, int y)							{ m_x = x; m_y = y; }
	int getX() const								{ return m_x; }
	int getY() const								{ return m_y; }
	int getComponent(int i) const							{ if (i == 1) return m_x; return m_y; }
	IntVector2D &operator+=(const IntVector2D v)					{ m_x += v.m_x; m_y += v.m_y; return *this; }
	IntVector2D &operator-=(const IntVector2D v)					{ m_x -= v.m_x; m_y -= v.m_y; return *this; }
	bool operator==(const IntVector2D v)						{ if (v.m_x == m_x && v.m_y == m_y) return true; return false; }
private:
	int m_x, m_y;
};

class IntLine2D
{
public:
	IntLine2D()									{ }
	IntLine2D(IntVector2D v, IntVector2D w)						{ m_v1 = v; m_v2 = w; }
	~IntLine2D()									{ }
	IntVector2D getVertex1() const							{ return m_v1; }
	IntVector2D getVertex2() const							{ return m_v2; }
	void setVertices(IntVector2D v, IntVector2D w)					{ m_v1 = v; m_v2 = w; }
private:
	IntVector2D m_v1, m_v2;
};

class Vector2DdPlus : public Vector2Dd
{
public:
	Vector2DdPlus(double x = 0,double y = 0,double v = 0) : Vector2D<double>(x,y)	{ m_value = v; }
	void setValue(double v)								{ m_value = v; }
	double getValue() const								{ return m_value; }
private:	
	double m_value;
};

class IntVector2DPlus : public IntVector2D
{
public:
	IntVector2DPlus(int x = 0, int y = 0, int v = 0) : IntVector2D(x,y)		{ m_value = v; }
	void getValue(int v)								{ m_value = v; }
	int getValue() const								{ return m_value; }
private:
	int m_value;
};

} // end namespace

#endif // VECTOR2D_H

