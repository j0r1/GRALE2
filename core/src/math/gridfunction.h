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
 * \file gridfunction.h
 */

#ifndef GRALE_GRIDFUNCTION_H

#define GRALE_GRIDFUNCTION_H

#include "graleconfig.h"
#include "real2dderivablefunction.h"
#include <vector>

namespace grale
{

/** Can be used to interpolate gridded data. */
class GRALE_IMPORTEXPORT GridFunction : public Real2DDerivableFunction
{
public:
	GridFunction(Real2DFunction &f ,Vector2D<double> bottomLeft, Vector2D<double> topRight, 
	             int numX = 512, int numY = 512, bool abs = true);
	GridFunction(const double *pValues, Vector2D<double> bottomLeft, Vector2D<double> topRight,
		     int numX, int numY, bool asPixels = false);
	~GridFunction();

	bool isValid()												{ return m_pValues != nullptr; }
	const double *getValues() const								{ return m_pValues; }
	int getWidth() const									{ return m_numX; }
	int getHeight() const									{ return m_numY; }
	Vector2D<double> getBottomLeft() const							{ return m_bottomLeft; }
	Vector2D<double> getTopRight() const							{ return m_topRight; }
	Vector2D<double> getPosition(int x, int y) const					{ return m_bottomLeft + Vector2D<double>(m_xStep*(double)x, m_yStep*(double)y); }
	double getValue(int x, int y) const							{ return m_pValues[x+y*m_numX]; }
	IntVector2D getPosition(Vector2D<double> v) const;

	void dump() const;
		
	double operator()(Vector2D<double> v) const;
	Vector2D<double> getGradient(Vector2D<double> v) const;
private:
	const double *GridFunctionConstructor(Real2DFunction &f ,Vector2D<double> bottomLeft, Vector2D<double> topRight,
										  int numX, int numY, bool abs);

	class PixelValue
	{
	public:
		PixelValue()									{ }
		PixelValue(int x, int y, double value)						{ m_x = x; m_y = y; m_value = value; }
		int getX() const								{ return m_x; }
		int getY() const								{ return m_y; }
		double getValue() const								{ return m_value; }
	private:
		int m_x, m_y;
		double m_value;
	};

	class MinMaxValue
	{
	public:
		MinMaxValue()									{ m_set = false; }
		void process(int x, double value)
		{
			if (!m_set)
			{
				m_set = true;
				m_minX = x;
				m_maxX = x;
				m_minValue = value;
				m_maxValue = value;
			}
			else
			{
				if (x < m_minX)
				{
					m_minX = x;
					m_minValue = value;
				}
				else if (x > m_maxX)
				{
					m_maxX = x;
					m_maxValue = value;
				}
			}
		}
		bool isSet() const								{ return m_set; }
		int getMinX() const								{ return m_minX; }
		int getMaxX() const								{ return m_maxX; }
		double getMinValue() const							{ return m_minValue; }
		double getMaxValue() const							{ return m_maxValue; }
	private:
		bool m_set;
		int m_minX, m_maxX;
		double m_minValue, m_maxValue;
	};
	
	Vector2D<double> m_topRight, m_bottomLeft;
	std::vector<double> m_values;
	const double *m_pValues;
	double m_xStep, m_yStep;
	int m_numX, m_numY;
	bool m_asPixels;
};

inline IntVector2D GridFunction::getPosition(Vector2D<double> v) const
{
	Vector2D<double> diff = v-m_bottomLeft;

	int xpos = (int)((diff.getX()/m_xStep)+0.5);
	int ypos = (int)((diff.getY()/m_yStep)+0.5);
	
	return IntVector2D(xpos, ypos);
}

} // end namespace

#endif // GRALE_GRIDFUNCTION_H
