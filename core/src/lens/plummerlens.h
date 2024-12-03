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
 * \file plummerlens.h
 */

#ifndef GRALE_PLUMMERLENS_H

#define GRALE_PLUMMERLENS_H

#include "graleconfig.h"
#include "symmetriclens.h"

namespace grale
{

/** Parameters for a lens with a Plummer density distribution. */
class GRALE_IMPORTEXPORT PlummerLensParams : public GravitationalLensParams
{
public:
	PlummerLensParams()								{ lensmass = 0; angwidth = 0; }

	/** Initializes the parameters.
	 *  Initializes the parameters.
	 *  \param mass Total mass of the lens.
	 *  \param angularwidth Angular width of the mass distribution.
	 */
	PlummerLensParams(double mass,double angularwidth)				{ lensmass = mass; angwidth = angularwidth; }

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

/** Describes a lens with a Plummer mass distribution. */
class GRALE_IMPORTEXPORT PlummerLens : public SymmetricLens
{
public:
	PlummerLens();
	~PlummerLens();

	std::unique_ptr<GravitationalLens> createUninitializedInstance() const override { return std::make_unique<PlummerLens>(); }

	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
	                           double *pPotentialValue) const;

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(double deflectionScale, double potentialScale, std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
	std::vector<CLFloatParamInfo> getCLAdjustableFloatingPointParameterInfo(double deflectionScale, double potentialScale) const override;
	std::unique_ptr<GravitationalLensParams> createLensParamFromCLFloatParams(double deflectionScale, double potentialScale, float *pFloatParams) const override;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
	bool getSurfaceMassDensityDerivative(double thetaLength, double &deriv) const override;
private:
	double mass,angularwidth,angularwidth2;
};

} // end namespace

#endif // GRALE_PLUMMERLENS_H

