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
 * \file squarelens.h
 */

#ifndef GRALE_SQUARELENS_H

#define GRALE_SQUARELENS_H

#include "graleconfig.h"
#include "gravitationallens.h"

namespace grale
{

/** Parameters for a lens with a square mass density distribution. */
class GRALE_IMPORTEXPORT SquareLensParams : public GravitationalLensParams
{
public:
	SquareLensParams()								{ lensmass = 0; angwidth = 0; }

	/** Initializes the parameters.
	 *  Initializes the parameters.
	 *  \param mass Total mass of the lens.
	 *  \param angularwidth Angular width of the mass distribution.
	 */
	SquareLensParams(double mass, double angularWidth)				{ lensmass = mass; angwidth = angularWidth; }

	/** Returns the currently set mass of the lens. */
	double getLensMass() const							{ return lensmass; }

	/** Returns the currently set angular width of the mass distribution. */	
	double getAngularWidth() const							{ return angwidth; }
	
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	double lensmass,angwidth;
};

/** Describes a lens with a square mass distribution. */
class GRALE_IMPORTEXPORT SquareLens : public GravitationalLens
{
public:
	SquareLens();
	~SquareLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(std::string &subRoutineName) const;
	int getCLSubLenses() const							{ return 1; }
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	static double deflectionFunction1(double x, double y)				{ return (0.5*y*LN(x*x+y*y)+x*ATAN(y/x)); }
	static double deflectionFunction2(double x, double y) 				{ return (0.5*x*LN(x*x+y*y)+y*ATAN(x/y)); }
	static double deriv11(double x, double y) 					{ return ATAN(y/x); }
	static double deriv22(double x, double y) 					{ return ATAN(x/y); }
	static double deriv12(double x, double y) 					{ return 0.5*LN(x*x+y*y); }
	static double potential(double x, double y) 					{ return 0.5*(x*x*ATAN(y/x)+y*y*ATAN(x/y)+x*y*(LN(x*x+y*y))); }

	double mass,angularwidth,angularwidth2;
	double angularscale,angularscale2;
	double relativeangularwidth2;
	double rawd2,dens;
};		      

} // end namespace

#endif // GRALE_SQUARELENS_H

