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
 * \file multiplegausslens.h
 */

#ifndef GRALE_MULTIPLEGAUSSLENS_H

#define GRALE_MULTIPLEGAUSSLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include "gausslensinfo.h"
#include <vector>

namespace grale
{

/** Parameters for a lens consisting of several Gauss distributions. */
class GRALE_IMPORTEXPORT MultipleGaussLensParams : public GravitationalLensParams
{
public:
	MultipleGaussLensParams()							{ }

	/** Initializes the parameters with those in the \c lensinfo list. */
	MultipleGaussLensParams(const std::vector<GaussLensInfo> &lensInfo)		{ m_lensInfo = lensInfo; }

	/** Adds a particular Gauss distribution. */
	void addLensInfo(const GaussLensInfo &inf) 					{ m_lensInfo.push_back(inf); }

	/** Returns the currently stored information. */
	const std::vector<GaussLensInfo> &getLensInfo() const				{ return m_lensInfo; }
	
	GravitationalLensParams *createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	std::vector<GaussLensInfo> m_lensInfo;
};

/** Describes a lens with a mass distribution which is the sum of several Gauss distributions. */
class GRALE_IMPORTEXPORT MultipleGaussLens : public GravitationalLens
{
public: 
	MultipleGaussLens();
	~MultipleGaussLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const ;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	double getDeriv11(Vector2D<double> theta) const;
	double getDeriv12(Vector2D<double> theta) const;
	double getDeriv22(Vector2D<double> theta) const;
	
	GaussLensInfo *lensinfo;
	double scalefactor;
	double totalmass;
	int numlenses;
};

} // end namespace

#endif // GRALE_MULTIPLEGAUSSLENS_H

