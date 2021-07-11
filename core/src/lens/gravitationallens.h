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
 * \file gravitationallens.h
 */

#ifndef GRALE_GRAVITATIONALLENS_H

#define GRALE_GRAVITATIONALLENS_H

#include "graleconfig.h"
#include "vector2d.h"
#include <serut/serializationinterface.h>
#include <string>
#include <memory>

namespace grale
{

/** Base class for gravitational lens parameters. */
class GRALE_IMPORTEXPORT GravitationalLensParams : public errut::ErrorBase
{
public:
	GravitationalLensParams()								{ }
	virtual ~GravitationalLensParams()							{ }

	/** Writes the paramters to a serut::SerializationInterface implementation. */
	virtual bool write(serut::SerializationInterface &si) const						{ setErrorString("Not implemented"); return false; }

	/** Reads parameters from a serut::SerializationInterface instance. */
	virtual bool read(serut::SerializationInterface &si)							{ setErrorString("Not implemented"); return false; }

	/** Creates a copy of the parameters. */
	virtual std::unique_ptr<GravitationalLensParams> createCopy() const = 0;
};

/** Base class for gravitational lens implementations. */
class GRALE_IMPORTEXPORT GravitationalLens : public errut::ErrorBase
{
public:
	/** Specific lens types. */
	enum LensType 
	{
		Gaussian, 		/**< A lens with a Gaussian density distribution. */
		MultiplePlummers, 	/**< A lens consisting of multiple Plummer distributions. */
		Plummer, 		/**< A lens with a Plummer density profile. */
		Pointmass, 		/**< A point mass lens. */
		SIS, 			/**< A Singular Isothermal Sphere (SIS). */
		NSIE, 			/**< A Non-Singular Isothermal Ellipse (NSIE). */
		NSIS, 			/**< A Non-Singular Isothermal Sphere (NSIS). */
		SIE, 			/**< A Singular Isothermal Ellipse (SIE). */
		Square,			/**< A lens with a square-shaped mass distribution. */
		MultipleSquares,	/**< A lens consisting of multiple square-shaped distributions. */
		MultipleGaussians,	/**< A lens consisting of multiple Gaussian density distributions. */
		MassSheet,		/**< A sheet of constant mass density. */
		Composite,		/**< A lens composed of other lenses. */
		MassDisk,		/**< A disk of constant mass density. */
		Profile,		/**< A circularly symmetric lens with a specific density profile. */
		PolynomialMassProfile,	/**< A circularly symmetric lens with a total mass profile composed of several polynomials. */
		MultipleWendland,	/**< A lens based on multiple Wendland W7,3 functions. */
		DeflectionGrid,         /**< A lens based on gridded deflection vector information. */
		NFW,			/**< A symmetric projected NFW profile. */
		EllipticNFW,		/**< An elliptical generalization of the symmetric NFW profile. */
		Sersic,			/**< A symmetric Sersic profile. */
		EllipticSersic,		/**< An elliptical generalization of the symmetric Sersic profile. */
		PIEMD, /** PIEMD (dPIE) as described in Eliasdottir (2014). */
		PIMD, /** Circularly symmetric version of PIEMD. */
		AlphaPot, /** Softened power law potential (similar to 'alphapot' in gravlens/lensmodel) */
		Harmonic, /** 2D harmonic mass density (cfr JPEG lensing paper) */
		PotentialGrid, /** Based on discrete set of projected potential values */
		CircularPieces, /** Exact lens potentials within circular regions, interpolated ones in between */
		MPContainer, /** Container for multiple lenses at different redshifts. */
		MaxLensType
	};
protected:
	/** Meant to be used by a specific lens implementation. */
	GravitationalLens(LensType t);
public:
	~GravitationalLens();

	/** Returns the type of the lens. */
	LensType getLensType() const								{ return m_lensType; }
	
	/** Initializes the lens with parameters \c params, setting its angular 
	 *  diameter distance to \c D_d.
	 */
	virtual bool init(double D_d, const GravitationalLensParams *pLensParams = 0);

	/** Calculates the result of the lens equation.
	 *  Calculates the result of the lens equation.
	 *  \param D_s Angular diameter distance between source and observer.
	 *  \param D_ds Angular diameter distance between source and lens.
	 *  \param theta theta-vector for which the lens equation should be calculated.
	 *  \param pBeta The result of the lens equation is stored here.
	 */
	bool traceTheta(double D_s, double D_ds, Vector2D<double> theta, Vector2D<double> *pBeta) const;

	/** Calculates the deflection angle for a given vector \c theta and stores the
	 *  result in \c alpha.
	 */
	virtual bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const = 0;

	/** Calculate derivatives of the deflection angle */
	virtual bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, 
			                       double &ayy, double &axy) const;

	/** Calculate inverse magnification from derivatives of deflection angle. */
	static double getInverseMagnification(double D_s, double D_ds, double axx, 
			                      double ayy, double axy) 				{ double f = D_ds/D_s; return (1.0-f*axx)*(1.0-f*ayy)-f*f*axy*axy; }

	/** Calculate shear size and angle from derivatives of the deflection angle. */
	static void getShearInfo(double D_s, double D_ds, double axx, double ayy, 
			         double axy, double *pShearAngle, double *pShearSize)		{ double f = D_ds/D_s; double g1 = 0.5*(axx-ayy); double g2 = axy; Vector2D<double> gVector(g1, g2); double g = gVector.getLength(); if (g2 >= 0) *pShearAngle = 0.5*ACOS(g1/g); else *pShearAngle = -0.5*ACOS(g1/g); *pShearSize = g*f; }

	/** Calculate shear size and angle from shear components. */
	static void getShearInfo(double gamma1, double gamma2, 
			         double *pShearAngle, double *pShearSize)			{ Vector2D<double> gVector(gamma1, gamma2); double g = gVector.getLength(); if (gamma2 >= 0) *pShearAngle = 0.5*ACOS(gamma1/g); else *pShearAngle = -0.5*ACOS(gamma1/g); *pShearSize = g; }


	/** Returns the inverse magnification factor for a specific direction.
	 *  Returns the inverse magnification factor for a specific direction.  
	 *  \param D_s Angular diameter distance between source and observer.
	 *  \param D_ds Angular diameter distance between source and lens.
	 *  \param theta The direction for which the factor should be calculated.
	 */
	virtual double getInverseMagnification(double D_s, double D_ds, Vector2D<double> theta) const;

	/** Returns the surface mass density for a specific direction. */
	virtual double getSurfaceMassDensity(Vector2D<double> theta) const = 0;
	
	/** Returns the D_d parameter used in the GravitationalLens::init function. */
	double getLensDistance() const								{ return m_Dd; }

	/** Changes the lens distance. */
	void setLensDistance(double D_d)							{ m_Dd = D_d; }

	bool getTimeDelay(double z_d, double D_s, double D_ds, Vector2D<double> theta, 
			  Vector2D<double> beta, double *pTimeDelay) const;
	virtual bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
			                   double *pPotentialValue) const;
	
	/** Returns a pointer to a copy of the lens parameters used in the GravitationalLens::init function. */
	const GravitationalLensParams *getLensParameters() const				{ return m_pParameters.get(); }

	/** Sets a distance scale to be used when estimating the derivatives of the 
	 *  function beta(theta) numerically.
	 */
	virtual void setDerivativeAngularDistanceScale(double distanceScale)			{ m_derivDistScale = ABS(distanceScale); }
	
	/** Writes the current lens to a serut::SerializationInterface instance. */
	bool write(serut::SerializationInterface &si) const;

	/** Reads a lens instance from a serut::SerializationInterface object and stores the lens in
	 *  \c lens; an error message is stored in \c errstr if the function is not
	 *  successful.
	 */
	static bool read(serut::SerializationInterface &si, std::unique_ptr<GravitationalLens> &pLens, std::string &errorString);

	/** Writes the current lens to a file. */
	bool save(const std::string &fileName) const;

	/** Loads a lens instance from a file and stores the lens in \c lens. an error 
	 *  message is stored in \c errstr if the function is not successful.
	 */
	static bool load(const std::string &fileName, std::unique_ptr<GravitationalLens> &pLens, std::string &errorString);

	/** Creates a copy of the current lens instance. */
	std::unique_ptr<GravitationalLens> createCopy() const;

	// TODO: experimental
	// OpenCL stuff
	
	virtual bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	virtual bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	virtual bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	
	std::string getCLLensProgram(std::string &subRoutineName, bool derivatives = true, bool potential = true) const;
	std::string getCLLensQuantitiesStructure(bool derivatives = true, bool potential = true) const;

	// like this:
	// LensQuantities subRoutineName(float2 coord, __global const int *pIntParams, __global const float *pFloatParams)
	//
	// typedef struct
	// {
	// 	float alphaX;
	// 	float alphaY;
	// 	float potential;
	// 	float axx;
	// 	float ayy;
	// 	float axy;
	// } LensQuantities;
	virtual std::string getCLProgram(std::string &subRoutineName, bool derivatives = true, bool potential = true) const;
protected:
	/** Specific lens implementations implement this function to process the parameters
	 *  specified in the GravitationalLens::init function.
	 */
	virtual bool processParameters(const GravitationalLensParams *params) = 0;
private:
	bool m_init;
	double m_Dd;
	LensType m_lensType;
	std::unique_ptr<GravitationalLensParams> m_pParameters;
	double m_derivDistScale;
};

} // end namespace

#endif // GRALE_GRAVITATIONALLENS_H

