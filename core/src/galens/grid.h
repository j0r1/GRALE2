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

/*
class GRALE_IMPORTEXPORT Grid : public errut::ErrorBase
{
public:
	Grid();
	Grid(const Grid &src);
	~Grid();

	Grid &operator=(const Grid &src);

	void addSquares(const std::list<GridSquare> &extraSquares);
	void getSquares(std::list<GridSquare> &gridSquares) const;

	bool buildUniform(Vector2D<double> center, double size, int subdiv);
	bool buildDensity(Real2DFunction &f, Vector2D<double> center, double size, double subdivfraction, 
			  int xdiv, int ydiv, bool abs, bool keeplarger);
	bool buildDensity(Real2DFunction &f, Vector2D<double> center, double size, int minsquares, int maxsquares, 
			  int xdiv, int ydiv, bool abs, bool keeplarger);
	bool buildGradient(Real2DDerivableFunction &f, Vector2D<double> center, double size, double subdivfraction,
			  int xdiv, int ydiv, bool keeplarger);
	bool buildDensityAndGradient(Real2DDerivableFunction &f, Vector2D<double> center, double size, 
			             double denssubdivfraction, double gradsubdivfraction,
				     int xdiv, int ydiv, bool abs, bool denskeeplarger,
				     bool gradkeeplarger);
	bool buildDensityAndGradient(Real2DDerivableFunction &f, Vector2D<double> center, double size, int minsquares,
			             int maxsquares, int xdiv, int ydiv, bool abs, bool keeplarger);
private:
	void copyFrom(const Grid &src);

	class RationalNumber
	{
	public:
		RationalNumber(int num, int denom)						{ a = num; b = denom; }
		~RationalNumber()								{ }
		int GetNumerator() const							{ return a; }
		int GetDenominator() const							{ return b; }
		double GetRealValue() const							{ return ((double)a)/((double)b); }
		bool operator==(const RationalNumber &r) const					{ if (a == r.a && b == r.b) return true; return false; }
	private:
		int a, b;
	};

	class RationalGridSquare
	{
	public:
		RationalGridSquare(RationalNumber xc, RationalNumber yc, RationalNumber s) : xcenter(xc),ycenter(yc),sz(s) 
												{ m_marked = false; }
		RationalNumber GetCenterX() const						{ return xcenter; }
		RationalNumber GetCenterY() const						{ return ycenter; }
		RationalNumber GetSize() const							{ return sz; }
		bool operator==(const RationalGridSquare &g) const				{ if (xcenter == g.xcenter && ycenter == g.ycenter && sz == g.sz) return true; return false; }

		bool isMarked() const								{ return m_marked; }
		void setMarker() 								{ m_marked = true; }
	private:
		RationalNumber xcenter,ycenter,sz;
		bool m_marked;
	};

	class GradientFunction : public Real2DFunction
	{
	public:
		GradientFunction(Real2DDerivableFunction &g) : f(g)				{ }
		~GradientFunction()								{ }
		double operator()(Vector2D<double> v) const					{ return f.getGradient(v).getLength(); }
	private:
		Real2DDerivableFunction &f;
	};

	std::list<RationalGridSquare> grid;
	Vector2D<double> gridcenter;
	double gridsize;

	std::list<GridSquare> m_extraSquares;
};
*/

} // end namespace

#endif // GRALE_GRID_H

