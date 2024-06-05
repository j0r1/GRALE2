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
 * \file multipleplummerlens.h
 */

#ifndef GRALE_MULTIPLEPLUMMERLENS_H

#define GRALE_MULTIPLEPLUMMERLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include "plummerlensinfo.h"
#include <vector>

namespace grale
{

/** Parameters for a lens consisting of several Plummer distributions. */
class GRALE_IMPORTEXPORT MultiplePlummerLensParams : public GravitationalLensParams
{
public:
	MultiplePlummerLensParams()							{ }

	/** Initializes the parameters with those in the \c lensinfo list. */
	MultiplePlummerLensParams(const std::vector<PlummerLensInfo> &lensInfo)		{ m_lensInfo = lensInfo; }

	/** Adds a particular Plummer distribution. */
	void addLensInfo(const PlummerLensInfo &inf) 					{ m_lensInfo.push_back(inf); }

	/** Returns the currently stored information. */
	const std::vector<PlummerLensInfo> &getLensInfo() const				{ return m_lensInfo; }
	
	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	std::vector<PlummerLensInfo> m_lensInfo;
};

/** Describes a lens with a mass distribution which is the sum of several Plummer distributions. */
class GRALE_IMPORTEXPORT MultiplePlummerLens : public GravitationalLens
{
public: 
	MultiplePlummerLens();
	~MultiplePlummerLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	bool getAlphaVectorSecondDerivatives(Vector2D<double> theta, double &axxx, double &ayyy, double &axxy, double &ayyx) const override;

	bool getSuggestedScales(double *pDeflectionScale, double *pPotentialScale) const;
	bool getCLParameterCounts(int *pNumIntParams, int *pNumFloatParams) const;
	bool getCLParameters(double deflectionScale, double potentialScale, int *pIntParams, float *pFloatParams) const;
	std::string getCLProgram(std::string &subRoutineName, bool derivatives = true, bool potential = true) const override;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double getDerivAlphaTheta(int i,int j,Vector2D<double> theta) const;
	
	std::vector<PlummerLensInfo> lensinfo;
	double scalefactor;
	double totalmass;
	int numlenses;
};

} // end namespace

#endif // GRALE_MULTIPLEPLUMMERLENS_H

