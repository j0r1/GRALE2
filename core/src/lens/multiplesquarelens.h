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
 * \file multiplesquarelens.h
 */

#ifndef GRALE_MULTIPLESQUARELENS_H

#define GRALE_MULTIPLESQUARELENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include "squarelensinfo.h"
#include <vector>

namespace grale
{

/** Parameters for a lens consisting of several square-shaped distributions. */
class GRALE_IMPORTEXPORT MultipleSquareLensParams : public GravitationalLensParams
{
public:
	MultipleSquareLensParams()							{ }

	/** Initializes the parameters with those in the \c lensinfo list. */
	MultipleSquareLensParams(const std::vector<SquareLensInfo> &lensInfo)		{ m_lensInfo = lensInfo; }

	/** Adds a particular square shaped distribution. */
	void addLensInfo(SquareLensInfo inf) 						{ m_lensInfo.push_back(inf); }

	/** Returns the currently stored information. */
	const std::vector<SquareLensInfo> &getLensInfo() const				{ return m_lensInfo; }
	
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	std::vector<SquareLensInfo> m_lensInfo;
};

/** Describes a lens with a mass distribution which is the sum of several square-shaped distributions. */
class GRALE_IMPORTEXPORT MultipleSquareLens : public GravitationalLens
{
public: 
	MultipleSquareLens();
	~MultipleSquareLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const; 
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	static double deflectionFunction1(double x, double y)				{ return (0.5*y*LN(x*x+y*y)+x*ATAN(y/x)); }
	static double deflectionFunction2(double x, double y) 				{ return (0.5*x*LN(x*x+y*y)+y*ATAN(x/y)); }
	static double deriv11(double x, double y) 					{ return ATAN(y/x); }
	static double deriv22(double x, double y) 					{ return ATAN(x/y); }
	static double deriv12(double x, double y) 					{ return 0.5*LN(x*x+y*y); }
	static double potential(double x, double y)		 			{ return 0.5*(x*x*ATAN(y/x)+y*y*ATAN(x/y)+x*y*(LN(x*x+y*y))); }

	SquareLensInfo *lensinfo;
	double scalefactor;
	double totalmass;
	int numlenses;
};

} // end namespace

#endif // GRALE_MULTIPLESQUARELENS_H

