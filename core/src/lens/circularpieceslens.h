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
 * \file circularpieceslens.h
 */

#ifndef GRALE_CIRCULARPIECESLENS_H

#define GRALE_CIRCULARPIECESLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include <vector>
#include <memory>

namespace grale
{

class GRALE_IMPORTEXPORT CircularPieceInfo
{
public:
	CircularPieceInfo(const std::shared_ptr<GravitationalLens> &lens,
			                 double startRadius, double endRadius,
							 double potentialScale, double potentialOffset)
		: m_lens(lens), 
		  m_startRadius(startRadius),
		  m_endRadius(endRadius),
		  m_potentialScale(potentialScale),
		  m_potentialOffset(potentialOffset)
	{
	}

	CircularPieceInfo(const CircularPieceInfo &c)
		: m_lens(c.m_lens), 
		  m_startRadius(c.m_startRadius),
		  m_endRadius(c.m_endRadius),
		  m_potentialScale(c.m_potentialScale),
		  m_potentialOffset(c.m_potentialOffset)
	{
	}

	const std::shared_ptr<GravitationalLens> &getLens() const								{ return m_lens; }
	double getStartRadius() const															{ return m_startRadius; }
	double getEndRadius() const																{ return m_endRadius; }
	double getPotentialScale() const														{ return m_potentialScale; }
	double getPotentialOffset() const														{ return m_potentialOffset; }
private:
	const std::shared_ptr<GravitationalLens> m_lens;
	double m_startRadius, m_endRadius;
	double m_potentialScale, m_potentialOffset;
};

class GRALE_IMPORTEXPORT CircularPiecesLensParams : public GravitationalLensParams
{
public:
	CircularPiecesLensParams()																{ }
	CircularPiecesLensParams(const std::vector<CircularPieceInfo> &pieces);

	const std::vector<CircularPieceInfo> &getPiecesInfo() const								{ return m_pieces; }

	GravitationalLensParams *createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	std::vector<CircularPieceInfo> m_pieces;
};

class GRALE_IMPORTEXPORT CircularPiecesLens : public GravitationalLens
{
public:
	CircularPiecesLens();
	~CircularPiecesLens();
	bool getAlphaVector(Vector2D<double> theta,Vector2D<double> *pAlpha) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	void getLenses(Vector2D<double> theta, int &idx1, int &idx2, double &f, double &df) const;

	std::vector<std::shared_ptr<GravitationalLens>> m_lenses;
	std::vector<double> m_startRadius, m_endRadius;
	std::vector<double> m_scale, m_potentialOffset;
};

} // end namespace

#endif // GRALE_CIRCULARPIECESLENS_H

