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
#include "vector2d.h"
#include "grid.h"
//#include "rectangularintegrator.h"

#include "debugnew.h"

namespace grale
{

	/*
	
Grid::Grid()
{
}

Grid::Grid(const Grid &src)
{
	copyFrom(src);
}

Grid::~Grid()
{
}

Grid &Grid::operator=(const Grid &src)
{
	copyFrom(src);
	return *this;
}

void Grid::addSquares(const std::list<GridSquare> &squares)
{
	for (auto it = squares.begin() ; it != squares.end() ; it++)
		m_extraSquares.push_back(*it);
}

void Grid::getSquares(std::list<GridSquare> &gridsquares) const
{
	gridsquares.clear();

	double leftx = gridcenter.getX() - gridsize/2.0;
	double bottomy = gridcenter.getY() - gridsize/2.0;
	
	for (auto it = grid.begin() ; it != grid.end() ; it++)
	{
		RationalGridSquare gs = (*it);
		RationalNumber centerx = gs.GetCenterX();
		RationalNumber centery = gs.GetCenterY();
		RationalNumber size = gs.GetSize();
		double cx = leftx + centerx.GetRealValue()*gridsize;
		double cy = bottomy + centery.GetRealValue()*gridsize;
		double sz = size.GetRealValue()*gridsize;
		gridsquares.push_back(GridSquare(Vector2D<double>(cx,cy),sz));
	}

	for (auto it = m_extraSquares.begin() ; it != m_extraSquares.end() ; it++)
		gridsquares.push_back(*it);
}

bool Grid::buildUniform(Vector2D<double> center, double size, int subdiv)
{
	if (subdiv < 1)
	{
		setErrorString("Illegal subdivision number");
		return false;
	}

	if (size <= 0)
	{
		setErrorString("Grid size must be positive");
		return false;
	}

	int denom = subdiv*2;
	
	grid.clear();
	gridcenter = center;
	gridsize = size;
	
	for (int i = 1 ; i < denom ; i += 2)
	{
		for (int j = 1 ; j < denom ; j += 2)
		{
			grid.push_back(RationalGridSquare(RationalNumber(i,denom),RationalNumber(j,denom),RationalNumber(1,subdiv)));
		}
	}

	return true;
}

bool Grid::buildDensity(Real2DFunction &f, Vector2D<double> center, double size, double subdivfraction,
			  int xdiv, int ydiv, bool abs, bool keeplarger)
{
	if (xdiv < 1)
	{
		setErrorString("Illegal integration subdivision size (x direction)");
		return false;
	}

	if (ydiv < 1)
	{
		setErrorString("Illegal integration subdivision size (y direction)");
		return false;
	}
	
	if (size <= 0)
	{
		setErrorString("Grid size must be positive");
		return false;
	}

	if (subdivfraction <= 0)
	{
		setErrorString("Invalid subdivision fraction value");
		return false;
	}

	grid.clear();
	gridcenter = center;
	gridsize = size;
	grid.push_back(RationalGridSquare(RationalNumber(1,2),RationalNumber(1,2),RationalNumber(1,1)));

	RectangularIntegrator rint(center-Vector2D<double>(size/2.0,size/2.0),center+Vector2D<double>(size/2.0,size/2.0),xdiv,ydiv);
	double total = rint.integrate(f,abs);
	double xintsize = size/(double)xdiv;
	double yintsize = size/(double)ydiv;
	double leftx = gridcenter.getX() - gridsize/2.0;
	double bottomy = gridcenter.getY() - gridsize/2.0;
	
	bool done = false;
	while (!done)
	{
		std::list<RationalGridSquare> newgrid;
		std::list<RationalGridSquare>::const_iterator it;
		bool gotsubdivide = false;
		
		for (it = grid.begin() ; it != grid.end() ; it++)
		{
			RationalGridSquare gs = (*it);
			RationalNumber centerx = gs.GetCenterX();
			RationalNumber centery = gs.GetCenterY();
			RationalNumber size = gs.GetSize();
			double cx = leftx + centerx.GetRealValue()*gridsize;
			double cy = bottomy + centery.GetRealValue()*gridsize;
			double sz = size.GetRealValue()*gridsize;
		
			if (gs.isMarked())
				newgrid.push_back(gs);
			else
			{
				int nx = (int)(sz/xintsize);
				int ny = (int)(sz/yintsize);
	
				if (nx < 1)
					nx = 1;
				if (ny < 1)
					ny = 1;
				RectangularIntegrator ri(Vector2D<double>(cx-sz/2.0,cy-sz/2.0),Vector2D<double>(cx+sz/2.0,cy+sz/2.0),nx,ny);
	
				double squareval = ri.integrate(f,abs);
				
				if (squareval > subdivfraction)
				{
					int xn = centerx.GetNumerator(); 
					int xd = centerx.GetDenominator();
					int yn = centery.GetNumerator();
					int yd = centery.GetDenominator();
					int s = size.GetDenominator();
					
					if (keeplarger)
					{
						gs.setMarker();
						newgrid.push_back(gs);
					}
					newgrid.push_back(RationalGridSquare(RationalNumber(xn*2-1,xd*2),RationalNumber(yn*2-1,yd*2),RationalNumber(1,2*s)));
					newgrid.push_back(RationalGridSquare(RationalNumber(xn*2-1,xd*2),RationalNumber(yn*2+1,yd*2),RationalNumber(1,2*s)));
					newgrid.push_back(RationalGridSquare(RationalNumber(xn*2+1,xd*2),RationalNumber(yn*2+1,yd*2),RationalNumber(1,2*s)));
					newgrid.push_back(RationalGridSquare(RationalNumber(xn*2+1,xd*2),RationalNumber(yn*2-1,yd*2),RationalNumber(1,2*s)));
	
					gotsubdivide = true;
				}
				else
				{
					gs.setMarker(); // don't need to recalculate on a next iteration
					newgrid.push_back(gs);
				}
			}
		}
		
		grid = newgrid;
		if (!gotsubdivide)
			done = true;
	}

	return true;
}

bool Grid::buildGradient(Real2DDerivableFunction &f, Vector2D<double> center, double size, double subdivfraction,
			  int xdiv, int ydiv, bool keeplarger)
{
	if (size <= 0)
	{
		setErrorString("Grid size must be positive");
		return false;
	}

	GradientFunction g(f);
	return buildDensity(g, center, size, subdivfraction, xdiv, ydiv, false, keeplarger);
}

bool Grid::buildDensityAndGradient(Real2DDerivableFunction &f, Vector2D<double> center, double size, 
			             double denssubdivfraction, double gradsubdivfraction,
				     int xdiv, int ydiv, bool abs, bool denskeeplarger,
				     bool gradkeeplarger)
{
	if (!buildDensity(f, center, size, denssubdivfraction, xdiv, ydiv, abs, denskeeplarger))
		return false;
	
	std::list<RationalGridSquare> denssquares = grid;
	std::list<RationalGridSquare>::const_iterator it;
	
	if (!buildGradient(f, center, size, gradsubdivfraction, xdiv, ydiv, gradkeeplarger))
		return false;
	
	// merge previous list
	
	std::list<RationalGridSquare> addsquares;
	
	for (it = denssquares.begin() ; it != denssquares.end() ; it++)
	{
		bool found = false;
		std::list<RationalGridSquare>::const_iterator it2 = grid.begin();

		while (!found && it2 != grid.end())
		{
			if ((*it2) == (*it))
				found = true;
			else
				it2++;
		}

		if (!found)
			addsquares.push_back(*it);
	}

	for (it = addsquares.begin() ; it != addsquares.end() ; it++)
		grid.push_back(*it);
	
	return true;
}

bool Grid::buildDensity(Real2DFunction &f, Vector2D<double> center, double size, int minsquares, int maxsquares, 
			int xdiv, int ydiv, bool abs, bool keeplarger)
{
	if (minsquares >= maxsquares)
	{
		setErrorString("The minimum number of squares specified is larger than the maximum "
			       "number of squares");
		return false;
	}
	
	bool done = false;
	double subdivfraction = 0.2;
	double diffrac = 2.0;
	
	for (int i = 0 ; !done && i < 50 ; i++) // make sure we end it some time
	{
		diffrac /= 2.0;
		
		do
		{
			subdivfraction /= (1.0+diffrac);
			if (!buildDensity(f, center, size, subdivfraction, xdiv, ydiv, abs, keeplarger))
				return false;
		} while (grid.size() < minsquares);

		if (grid.size() <= maxsquares)
			done = true;
		else
		{
			diffrac /= 2.0;
			
			do
			{
				subdivfraction *= (1.0+diffrac);
				if (!buildDensity(f, center, size, subdivfraction, xdiv, ydiv, abs, keeplarger))
					return false;
			} while (grid.size() > maxsquares);

			if (grid.size() >= minsquares)
				done = true;
		}
	}

	if (!done)
	{
		grid.clear();
		setErrorString("A grid of the requested size could not be constructed");
		return false;
	}
	
	return true;
}

bool Grid::buildDensityAndGradient(Real2DDerivableFunction &f, Vector2D<double> center, double size, int minsquares, int maxsquares, 
			int xdiv, int ydiv, bool abs, bool keeplarger)
{
	if (minsquares >= maxsquares)
	{
		setErrorString("The minimum number of squares specified is larger than the maximum "
			       "number of squares");
		return false;
	}
	
	bool done = false;
	double subdivfraction = 0.2;
	double diffrac = 2.0;
	
	for (int i = 0 ; !done && i < 50 ; i++) // make sure we end it some time
	{
		diffrac /= 2.0;
		
		do
		{
			subdivfraction /= (1.0+diffrac);
			if (!buildDensityAndGradient(f, center, size, subdivfraction, subdivfraction, xdiv, ydiv, abs, keeplarger, keeplarger))
				return false;
		} while (grid.size() < minsquares);

		if (grid.size() <= maxsquares)
			done = true;
		else
		{
			diffrac /= 2.0;
			
			do
			{
				subdivfraction *= (1.0+diffrac);
				if (!buildDensityAndGradient(f, center, size, subdivfraction, subdivfraction, xdiv, ydiv, abs, keeplarger, keeplarger))
					return false;
			} while (grid.size() > maxsquares);

			if (grid.size() >= minsquares)
				done = true;
		}
	}

	if (!done)
	{
		grid.clear();
		setErrorString("A grid of the requested size could not be constructed");
		return false;
	}
	
	return true;
}

void Grid::copyFrom(const Grid &src)
{
	grid = src.grid;
	gridcenter = src.gridcenter;
	gridsize = src.gridsize;

	m_extraSquares = src.m_extraSquares;
}
*/

} // end namespace

