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
 * \file polynomialmassprofilelens.h
 */

#ifndef GRALE_POLYNOMIALMASSPROFILELENS_H

#define GRALE_POLYNOMIALMASSPROFILELENS_H

#include "graleconfig.h"
#include "symmetriclens.h"
#include <vector>

// Here, the profile refers to the profile of the total mass M(theta)

namespace grale
{

class GRALE_IMPORTEXPORT PolynomialPart
{
public:
	PolynomialPart();
	PolynomialPart(double xOffset, double yOffset, double xScale, double yScale, double xEnd, const std::vector<double> &coefficients);
	PolynomialPart(const PolynomialPart &polynomialPart);
	~PolynomialPart();
	double getXOffset() const									{ return m_xOffset; }
	double getYOffset() const									{ return m_yOffset; }
	double getXScale() const									{ return m_xScale; }
	double getYScale() const									{ return m_yScale; }
	double getEndPosition() const									{ return m_xEnd; }
	const std::vector<double> &getCoefficients() const						{ return m_coefficients; }
	const std::vector<double> &getOffsetCoefficients() const					{ return m_offsetCoefficients; }

	PolynomialPart &operator=(const PolynomialPart &p);
private:
	void calculateOffsetCoefficients();
	
	double m_xOffset;
	double m_yOffset;
	double m_xScale;
	double m_yScale;
	double m_xEnd;
	std::vector<double> m_coefficients;
	std::vector<double> m_offsetCoefficients;
};

class GRALE_IMPORTEXPORT PolynomialMassProfileLensParams : public GravitationalLensParams
{
public:
	PolynomialMassProfileLensParams();
	~PolynomialMassProfileLensParams();
	void addPolynomialPart(const PolynomialPart &p);
	const std::vector<PolynomialPart> &getPolynomialParts() const					{ return m_polynomialParts; }

	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
	GravitationalLensParams *createCopy() const;
private:
	std::vector<PolynomialPart> m_polynomialParts;
};

class GRALE_IMPORTEXPORT PolynomialZeroMassLensParams : public PolynomialMassProfileLensParams
{
public:
	PolynomialZeroMassLensParams(double D_d, double densityScale, double angularRadius, double zeroPoint);
	~PolynomialZeroMassLensParams();
};

class GRALE_IMPORTEXPORT PolynomialTimeDelayLensParams : public PolynomialMassProfileLensParams
{
public:
	PolynomialTimeDelayLensParams(double theta1, double theta2, double timeDiff, double zLens);
	~PolynomialTimeDelayLensParams();
};

/** Bla
 *
 * For the potential, we'll need to calculate something like:
 * \f[
 * 	\int_{x_1}^{x_2}\frac{\sum_{k=0}^M a_k x^k}{x+b} dx
 * \f]
 *
 * To calculate this, the following relation can be helpful
 * \f[
 * 	\sum_{k=0}^M a_k x^k = \sum_{k=0}^M f_k (x+b)^k
 * \f]
 * in which
 * \f[
 * 	f_k = \sum_{l=k}^M \left(\begin{array}{c} l \\ k \end{array}\right) a_l (-b)^{l-k}
 * \f]
 * This can be proved by noting that
 * \f[
 * 	x^N = [(x+b)-b]^N = \sum_{k=0}^N \left(\begin{array}{c} N \\ k \end{array}\right) (x+b)^k (-b)^{N-k}
 * \f]
 * so that
 * \f[
 *	\sum_{l=0}^M a_l x^l = \sum_{l=0}^M \sum_{k=0}^l \left(\begin{array}{c} l \\ k \end{array}\right) a_l (x+b)^k (-b)^{l-k}
 * \f]
 
 */

class GRALE_IMPORTEXPORT PolynomialMassProfileLens : public SymmetricLens
{
public:
	PolynomialMassProfileLens();
	~PolynomialMassProfileLens();
	
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
	                           double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
	double getMassInside(double thetaLength) const;
	double getProfileSurfaceMassDensity(double thetaLength) const;
private:
	std::vector<PolynomialPart> m_polynomialParts;
	std::vector<double> m_potentialOffsets;
	double m_totalMass;
	double m_totalPotential;
	double m_endPosition;
};

} // end namespace

#endif // GRALE_POLYNOMIALMASSPROFILELENS_H

