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
 * \file polygon2d.h
 */

#ifndef GRALE_POLYGON2D_H

#define GRALE_POLYGON2D_H

#include "graleconfig.h"
#include "vector2d.h"
#include "polygon2d.h"
#include "rectangle2d.h"
#include "line2d.h"
#include "constants.h"
#include "Wm5ConvexHull2.h"
#include "Wm5ContMinBox2.h"
#include <iterator>
#include <vector>

namespace grale
{

/** This class can be used to represent a polygon. */
template<class T>
class Polygon2D
{
public:
	Polygon2D()										{ m_numCoords = 0; m_isConvexHull = false; }
	Polygon2D(const Polygon2D<T> &p)							{ m_numCoords = 0; m_isConvexHull = false; copyFieldsFrom(p); }
	~Polygon2D()										{ }

	void init(const Vector2D<T> *pPoints, int numPoints, bool calcConvexHull);
	void init(const std::vector<Vector2D<T> > &points, bool calcConvexHull);
	bool init(Polygon2D<T> &hull, int numPoints);
	bool isInside(Vector2D<T> p) const;
	T getDistanceSquared(Vector2D<T> p) const;
	T getArea() const;
	bool isConvexHull() const								{ return m_isConvexHull; }
	const Vector2D<T> *getPoints() const							{ return &(m_xyCoords[0]); }
	int getNumberOfPoints() const								{ return m_numCoords; }
	void scale(T scaleFactor);
	
	bool getMinimalAreaRectangle(Rectangle2D<T> &r) const;
	
	void operator=(const Polygon2D<T> &p)							{ copyFieldsFrom(p); }
private:
	void copyFieldsFrom(const Polygon2D<T> &p);

	template<class Iterator> void calculateHull(Iterator startIt, int numPoints);
	static void calcRectangle(const Line2D<T> supportLines[4], Rectangle2D<T> &rect);
	inline T calcDist2(Vector2D<T> p1, Vector2D<T> p2, Vector2D<T> P) const;
	
	std::vector<Vector2D<T> > m_xyCoords;
	int m_numCoords;
	bool m_isConvexHull;
};

template<class T>
inline bool Polygon2D<T>::isInside(Vector2D<T> p) const
{
	T x = p.getX();
	T y = p.getY();
	int intersections = 0;
	
	for (int i = 0 ; i < m_numCoords ; i++)
	{
		T y1 = m_xyCoords[i].getY();
		T y2 = m_xyCoords[i+1].getY();
		T Y1 = MIN(m_xyCoords[i].getY(), m_xyCoords[i+1].getY());
		T Y2 = MAX(m_xyCoords[i].getY(), m_xyCoords[i+1].getY());

		if (Y1 < y && y <= Y2)
		{
			T x1 = m_xyCoords[i].getX();
			T x2 = m_xyCoords[i+1].getX();

			if (x <= MAX(x1, x2))
			{
				if (x1 == x2)
					intersections++;
				else
				{
					T x0 = ((y - y1)*(x2 - x1))/(y2 - y1) + x1;

					if (x <= x0)
						intersections++;
				}
			}
		}
	}
	return (intersections&1)?true:false;
}

template<class T>
inline T Polygon2D<T>::getDistanceSquared(Vector2D<T> p) const
{
	if (m_numCoords == 0)
		return 0;
	if (m_numCoords == 1)
	{
		Vector2D<T> p2 = m_xyCoords[0];
		p2 -= p;
		return p2.getLengthSquared();
	}
	
	T minDist2 = calcDist2(m_xyCoords[0], m_xyCoords[1], p);
	
	for (int i = 1 ; i < m_numCoords ; i++)
	{
		T d2 = calcDist2(m_xyCoords[i], m_xyCoords[i+1], p);

		if (d2 < minDist2)
			minDist2 = d2;
	}
	return minDist2;
}
		
template<class T>
inline T Polygon2D<T>::calcDist2(Vector2D<T> p1, Vector2D<T> p2, Vector2D<T> P) const
{
	Vector2D<T> Q;
	T dx = p2.getX() - p1.getX();
	T dy = p2.getY() - p1.getY();
	T denom = dx*dx + dy*dy;

	if (denom == 0)
		Q = p1;
	else
	{
		T l = ((P.getX() - p1.getX())*dx + (P.getY() - p1.getY())*dy)/denom;
		
		if (l <= 0)
			Q = p1;
		else if (l >= 1)
			Q = p2;
		else
			Q = Vector2D<T>(p1.getX() + l*dx, p1.getY() + l*dy);
	}
	
	Q -= P;
	return Q.getLengthSquared();
}
	
template<class T>
inline void Polygon2D<T>::copyFieldsFrom(const Polygon2D<T> &p)
{
	m_numCoords = p.m_numCoords;
	m_isConvexHull = p.m_isConvexHull;
	m_xyCoords = p.m_xyCoords;
}

template<class T>
inline T Polygon2D<T>::getArea() const
{
	T area = 0;

	for (int i = 0 ; i < m_numCoords ; i++)
	{
		T areaPart = m_xyCoords[i].getX()*m_xyCoords[i+1].getY() - m_xyCoords[i+1].getX()*m_xyCoords[i].getY();

		area += areaPart;
	}

	area /= (T)2;

	return area;
}

template<class T>
inline void Polygon2D<T>::scale(T scaleFactor)
{
#if 0
	Vector2D<T> center(0, 0);
	for (int i = 0 ; i < m_numCoords ; i++)
		center += m_pXYCoords[i];
	center /= (T)m_numCoords;
#else
	
	T area = 0;
	T cx = 0;
	T cy = 0;

	for (int i = 0 ; i < m_numCoords ; i++)
	{
		T areaPart = m_xyCoords[i].getX()*m_xyCoords[i+1].getY() - m_xyCoords[i+1].getX()*m_xyCoords[i].getY();

		area += areaPart;
		cx += (m_xyCoords[i].getX() + m_xyCoords[i+1].getX())*areaPart;
		cy += (m_xyCoords[i].getY() + m_xyCoords[i+1].getY())*areaPart;

	}

	area /= (T)2;
	cx /= area*(T)6;
	cy /= area*(T)6;

	Vector2D<T> center(cx, cy);

	for (int i = 0 ; i < m_numCoords ; i++)
	{
		Vector2D<T> diff = m_xyCoords[i] - center;
		diff *= scaleFactor;
		m_xyCoords[i] = center + diff;
	}
#endif
}

template<class T>
void Polygon2D<T>::init(const Vector2D<T> *pPoints, int numPoints, bool calcConvexHull)
{
	if (numPoints == 0)
	{
		m_numCoords = 0;
		m_isConvexHull = false;
		m_xyCoords.resize(0);
		return;
	}
	
	if (calcConvexHull)
	{
		if (numPoints <= 3)
			init(pPoints, numPoints, false);
		else
			calculateHull(pPoints, numPoints);
		m_isConvexHull = true;
		return;
	}

	m_isConvexHull = false;

	m_xyCoords.resize(numPoints+1);
	m_numCoords = numPoints;
	
	const Vector2D<T> *pPtr;
	int i;

	pPtr = pPoints;
	
	m_xyCoords[m_numCoords] = *pPtr;
	for (i = 0 ;  i < m_numCoords ; pPtr++, i++)
		m_xyCoords[i] = *pPtr;
}

template<class T>
void Polygon2D<T>::init(const std::vector<Vector2D<T> > &points, bool calcConvexHull)
{
	if (points.size() == 0)
	{
		m_numCoords = 0;
		m_isConvexHull = false;
		m_xyCoords.resize(0);
		return;
	}
	
	if (calcConvexHull)
	{
		if (points.size() <= 3)
			init(points, false);
		else
			calculateHull(points.begin(), points.size());
		m_isConvexHull = true;
		return;
	}

	m_isConvexHull = false;

	m_xyCoords.resize(points.size() + 1);
	m_numCoords = points.size();
	
	typename std::vector< Vector2D<T> >::const_iterator it;
	int i;

	it = points.begin();
	
	m_xyCoords[m_numCoords] = (*it);
	for (i = 0 ;  i < m_numCoords ; it++, i++)
		m_xyCoords[i] = (*it);
}

template<class T>
template<class Iterator>	
void Polygon2D<T>::calculateHull(Iterator startIt, int numPoints)
{
	std::vector<Wm5::Vector2<T> > points(numPoints);

	Iterator it = startIt;

	for (int i = 0 ; i < numPoints ; i++, it++)
		points[i] = Wm5::Vector2<T>((*it).getX(), (*it).getY());

	m_xyCoords.resize(0);

	Wm5::ConvexHull2<T> hull(numPoints, &(points[0]), 0, false, Wm5::Query::QT_INTEGER);

	int numIndices = hull.GetNumSimplices();
	const int *pIndices = hull.GetIndices();

	m_numCoords = numIndices;
	m_xyCoords.resize(m_numCoords+1);

	m_xyCoords[m_numCoords] = Vector2D<T>(points[pIndices[0]].X(), points[pIndices[0]].Y());
	for (int i = 0 ; i < m_numCoords ; i++)
		m_xyCoords[i] = Vector2D<T>(points[pIndices[i]].X(), points[pIndices[i]].Y());
}
	
template<class T>
bool Polygon2D<T>::getMinimalAreaRectangle(Rectangle2D<T> &r) const
{
	if (!m_isConvexHull || m_numCoords < 3)
		return false;

	std::vector<Wm5::Vector2<T> > points(m_numCoords);

	for (int i = 0 ; i < m_numCoords ; i++)
		points[i] = Wm5::Vector2<T>(m_xyCoords[i].getX(), m_xyCoords[i].getY());

	Wm5::MinBox2<T> box(m_numCoords, &(points[0]), 0, Wm5::Query::QT_INTEGER, true);
	Wm5::Vector2<T> rectPoints[4];

	Wm5::Box2<T> rectBox = box;
	rectBox.ComputeVertices(rectPoints);
	
	Vector2D<T> graleRectPoints[4];

	graleRectPoints[0] = Vector2D<T>(rectPoints[0].X(), rectPoints[0].Y());
	graleRectPoints[1] = Vector2D<T>(rectPoints[1].X(), rectPoints[1].Y());
	graleRectPoints[2] = Vector2D<T>(rectPoints[2].X(), rectPoints[2].Y());
	graleRectPoints[3] = Vector2D<T>(rectPoints[3].X(), rectPoints[3].Y());
	r.init(graleRectPoints);
	return true;
}

typedef Polygon2D<double> Polygon2Dd;
typedef Polygon2D<float> Polygon2Df;

} // end namespace 

#endif // GRALE_POLYGON2D_H

