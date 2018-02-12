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
 * \file triangle2d.h
 */

#ifndef GRALE_TRIANGLE2D_H

#define GRALE_TRIANGLE2D_H

#include "graleconfig.h"
#include "line2d.h"
#include "Wm5IntrTriangle2Triangle2.h"
#include <algorithm>
#include <vector>

namespace grale
{

/** This class can be used to represent a triangle. */
template<class T>
class Triangle2D
{
public:
	Triangle2D()											{ }
	Triangle2D(Vector2D<T> point1, Vector2D<T> point2, Vector2D<T> point3)				{ m_points[0] = point1; m_points[1] = point2; m_points[2] = point3; m_points[3] = point1; }
	Triangle2D(Vector2D<T> p[3])									{ m_points[0] = p[0]; m_points[1] = p[1]; m_points[2] = p[2]; m_points[3] = p[0]; }
	T getArea() const;
	T getOverlapArea(const Triangle2D<T> &t) const;
	bool isInside(Vector2D<T> point) const;
	const Vector2D<T> *getPoints() const								{ return m_points; }
	Vector2D<T> getCentroid() const;
protected:
	Vector2D<T> m_points[4]; // last one is same as first one
};

typedef Triangle2D<double> Triangle2Dd;
typedef Triangle2D<float> Triangle2Df;

template<class T>
inline T Triangle2D<T>::getArea() const
{
	// Easy to calculate using cross product
	
	Vector2D<T> vec1 = m_points[1] - m_points[0];
	Vector2D<T> vec2 = m_points[2] - m_points[0];

	return ((T)0.5)*ABS(vec1.getX()*vec2.getY() - vec1.getY()*vec2.getX());
}

template<class T>
bool Triangle2D<T>::isInside(Vector2D<T> point) const
{
	int posCount = 0;
	int negCount = 0;
	
	for (int i = 0 ; i < 3 ; i++)
	{
		Vector2D<T> v = m_points[i+1] - m_points[i];
		Vector2D<T> w = point - m_points[i];
		T value = v.getX()*w.getY() - v.getY()*w.getX();
		if (value > 0)
			posCount++;
		else if (value < 0)
			negCount++;
		else
		{
			posCount++;
			negCount++;
		}
	}
	if (posCount == 3 || negCount == 3)
		return true;
	return false;
}

template<class T>
inline T Triangle2D<T>::getOverlapArea(const Triangle2D<T> &t) const
{
	Wm5::Triangle2<T> t1;
	Wm5::Triangle2<T> t2;

	Vector2D<T> A = m_points[1]-m_points[0];
	Vector2D<T> B = m_points[2]-m_points[0];

	if (A.getX()*B.getY() - A.getY()*B.getX() > 0)
	{
		t1 = Wm5::Triangle2<T>(Wm5::Vector2<T>(m_points[0].getX(), m_points[0].getY()),
			         Wm5::Vector2<T>(m_points[1].getX(), m_points[1].getY()),
			         Wm5::Vector2<T>(m_points[2].getX(), m_points[2].getY()));
	}
	else
	{
		t1 = Wm5::Triangle2<T>(Wm5::Vector2<T>(m_points[1].getX(), m_points[1].getY()),
			         Wm5::Vector2<T>(m_points[0].getX(), m_points[0].getY()),
			         Wm5::Vector2<T>(m_points[2].getX(), m_points[2].getY()));
	}
	
	A = t.m_points[1]-t.m_points[0];
	B = t.m_points[2]-t.m_points[0];

	if (A.getX()*B.getY() - A.getY()*B.getX() > 0)
	{
		t2 = Wm5::Triangle2<T>(Wm5::Vector2<T>(t.m_points[0].getX(), t.m_points[0].getY()),
			         Wm5::Vector2<T>(t.m_points[1].getX(), t.m_points[1].getY()),
			         Wm5::Vector2<T>(t.m_points[2].getX(), t.m_points[2].getY()));
	}
	else
	{
		t2 = Wm5::Triangle2<T>(Wm5::Vector2<T>(t.m_points[1].getX(), t.m_points[1].getY()),
			         Wm5::Vector2<T>(t.m_points[0].getX(), t.m_points[0].getY()),
			         Wm5::Vector2<T>(t.m_points[2].getX(), t.m_points[2].getY()));
	}

	Wm5::IntrTriangle2Triangle2<T> intersection(t1, t2);

	//for (int i = 0 ; i < 4 ; i++)
		//std::cout << m_points[i].getX() << " " << m_points[i].getY() << std::endl;

	//std::cout << std::endl;
	//std::cout << std::endl;

	//for (int i = 0 ; i < 4 ; i++)
		//std::cout << t.m_points[i].getX() << " " << t.m_points[i].getY() << std::endl;

	//std::cout << std::endl;
	//std::cout << std::endl;

	if (!intersection.Test())
		return 0;

	if (!intersection.Find())
		return 0;

	int numIntersectionPoints = intersection.GetQuantity();

	//for (int i = 0 ; i < numIntersectionPoints ; i++)
		//std::cout << intersection.GetPoint(i).X() << " " << intersection.GetPoint(i).Y() << std::endl;
	//std::cout << intersection.GetPoint(0).X() << " " << intersection.GetPoint(0).Y() << std::endl;

	if (numIntersectionPoints < 3)
		return 0;

	int num = numIntersectionPoints - 2;
	T totalArea = 0;
	Vector2D<T> v1(intersection.GetPoint(0).X(), intersection.GetPoint(0).Y());

	for (int i = 0 ; i < num ; i++)
	{
		Vector2D<T> v2(intersection.GetPoint(i+1).X(),intersection.GetPoint(i+1).Y());
		Vector2D<T> v3(intersection.GetPoint(i+2).X(),intersection.GetPoint(i+2).Y());

		Triangle2D<T> t2(v1, v2, v3);

		totalArea += t2.getArea();
	}
	
	return totalArea;
}

template<class T>
Vector2D<T> Triangle2D<T>::getCentroid() const
{
	return (m_points[0] + m_points[1] + m_points[2])/((T)3);
}

template<class T>
class Triangle2DPlus : public Triangle2D<T>
{
public:
	Triangle2DPlus(Vector2D<T> point1, Vector2D<T> point2, Vector2D<T> point3, T height1, T height2, T height3) : Triangle2D<T>(point1, point2, point3)
	{
		m_heights[0] = height1;
		m_heights[1] = height2;
		m_heights[2] = height3;
		m_heights[3] = height1;
		m_diff1 = Triangle2D<T>::m_points[1] - Triangle2D<T>::m_points[0];
		m_diff2 = Triangle2D<T>::m_points[2] - Triangle2D<T>::m_points[0];
		m_hdiff1 = m_heights[1] - m_heights[0];
		m_hdiff2 = m_heights[2] - m_heights[0];
		m_denom = m_diff1.getX()*m_diff2.getY() - m_diff1.getY()*m_diff2.getX();
	}

	bool isInside(Vector2D<T> point, T &height)
	{
		if (!Triangle2D<T>::isInside(point))
			return false;
		
		Vector2D<T> pdiff = point - Triangle2D<T>::m_points[0];
		T alpha = (m_diff2.getY()*pdiff.getX() - m_diff2.getX()*pdiff.getY())/m_denom;
		T beta = (-m_diff1.getY()*pdiff.getX() + m_diff1.getX()*pdiff.getY())/m_denom;

		height = m_heights[0] + alpha*m_hdiff1  + beta*m_hdiff2;

		return true;
	}

	const double *getHeights() const						{ return m_heights; }
private:
	T m_heights[4];
	T m_denom, m_hdiff1, m_hdiff2;
	Vector2D<T> m_diff1, m_diff2;
};

typedef Triangle2DPlus<double> Triangle2DPlusd;
typedef Triangle2DPlus<float> Triangle2DPlusf;

}

#endif // GRALE_TRIANGLE2D_H

