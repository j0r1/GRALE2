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

#include "graleconfig.h"
#include "gridfunction.h"
#include <iostream>

#include "debugnew.h"

using namespace std;

namespace grale
{

const double *GridFunction::GridFunctionConstructor(Real2DFunction &f, Vector2D<double> bottomLeft, Vector2D<double> topRight,
									  int numX, int numY, bool abs)
{
	if (numX < 10)
		m_numX = 10;
	else
		m_numX = numX;
	if (numY < 10)
		m_numY = 10;
	else
		m_numY = numY;

	int numTotal = m_numX*m_numY;
	m_values.resize(numTotal);

	m_xStep = (topRight.getX() - bottomLeft.getX())/((double)(m_numX-1));
	m_yStep = (topRight.getY() - bottomLeft.getY())/((double)(m_numY-1));
	m_topRight = topRight;
	m_bottomLeft = bottomLeft;
	
	double x2, y2;
	int x, y;
	int i = 0;

	if (!abs)
	{
		for (y = 0, y2 = m_bottomLeft.getY() ; y < m_numY ; y++, y2 += m_yStep)
		{
			for (x = 0, x2 = m_bottomLeft.getX() ; x < m_numX ; x++, x2 += m_xStep)
			{
				m_values[i] = f(Vector2D<double>(x2, y2));
				i++;
			}
		}
	}
	else
	{
		for (y = 0, y2 = m_bottomLeft.getY() ; y < m_numY ; y++, y2 += m_yStep)
		{
			for (x = 0, x2 = m_bottomLeft.getX() ; x < m_numX ; x++, x2 += m_xStep)
			{
				m_values[i] = ABS(f(Vector2D<double>(x2, y2)));
				i++;
			}
		}	
	}

	return &m_values[0];
}

GridFunction::GridFunction(Real2DFunction &f, Vector2D<double> bottomLeft, Vector2D<double> topRight,
                           int numX, int numY, bool abs) 
	: m_pValues(GridFunctionConstructor(f, bottomLeft, topRight, numX, numY, abs)),
	  m_asPixels(false)
{
}

GridFunction::GridFunction(const double *pValues, Vector2D<double> bottomLeft, Vector2D<double> topRight,
		           int numX, int numY, bool asPixels) 
	: m_pValues(pValues), m_asPixels(asPixels)
{
	m_numX = numX;
	m_numY = numY;

	m_pValues = (double *)pValues;

	if (!asPixels)
	{
		m_xStep = (topRight.getX() - bottomLeft.getX())/((double)(m_numX-1));
		m_yStep = (topRight.getY() - bottomLeft.getY())/((double)(m_numY-1));
	}
	else
	{
		m_xStep = (topRight.getX() - bottomLeft.getX())/((double)(m_numX));
		m_yStep = (topRight.getY() - bottomLeft.getY())/((double)(m_numY));
	}

	m_topRight = topRight;
	m_bottomLeft = bottomLeft;
}

GridFunction::~GridFunction()
{
}

double GridFunction::operator()(Vector2D<double> v) const
{
	double xDist = v.getX() - m_bottomLeft.getX();
	double yDist = v.getY() - m_bottomLeft.getY();
	
	double xFracDist = xDist/m_xStep;
	double yFracDist = yDist/m_yStep;
	
	int X = (int)(xFracDist);
	int Y = (int)(yFracDist);

	if (m_asPixels)
	{
		if (X < 0) X = 0;
		else if (X >= m_numX) X = m_numX-1;
		if (Y < 0) Y = 0;
		else if (Y >= m_numY) Y = m_numY-1;

		return m_pValues[X+Y*m_numX];
	}

	double t = xFracDist-(double)X;
	double u = yFracDist-(double)Y;

	if (xFracDist < 0) 
	{
		X--;
		t = 1.0-t;
	}
	if (yFracDist < 0)
	{
		Y--;
		u = 1.0-u;
	}

	int X0 = X;
	int X1 = X+1;
	int Y0 = Y;
	int Y1 = Y+1;

	if (X0 < 0) X0 = 0;
	if (X1 < 0) X1 = 0;
	if (X0 >= m_numX) X0 = m_numX-1;
	if (X1 >= m_numX) X1 = m_numX-1;

	if (Y0 < 0) Y0 = 0;
	if (Y1 < 0) Y1 = 0;
	if (Y0 >= m_numY) Y0 = m_numY-1;
	if (Y1 >= m_numY) Y1 = m_numY-1;
	
	double y1 = m_pValues[X0 + m_numX*Y0];
	double y2 = m_pValues[X1 + m_numX*Y0];
	double y3 = m_pValues[X1 + m_numX*Y1];
	double y4 = m_pValues[X0 + m_numX*Y1];

	double val = (1.0-t)*(1.0-u)*y1 + t*(1.0-u)*y2 + t*u*y3 + (1.0-t)*u*y4;
	return val;
}

Vector2D<double> GridFunction::getGradient(Vector2D<double> v) const
{
	double a, x, xDiff;
	double y, yDiff;

	xDiff = m_xStep/10.0;
	yDiff = m_yStep/10.0;

	a = operator()(v);
	x = operator()(v + Vector2D<double>(xDiff, 0.0));
	y = operator()(v + Vector2D<double>(0.0, yDiff));
	return Vector2D<double>((x-a)/xDiff, (y-a)/yDiff);
}

void GridFunction::dump() const
{
	int idx = 0;
	
	for (int y = 0 ; y < m_numY ; y++)
	{
		for (int x = 0 ; x < m_numX ; x++, idx++)
		{
			std::cout << x << " " << y << " " << m_pValues[idx] << std::endl;
		}
		std::cout << std::endl;
	}
}

} // end namespace

