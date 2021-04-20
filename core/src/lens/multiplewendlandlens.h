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

#ifndef GRALE_MULTIPLEWENDLANDLENS_H

#define GRALE_MULTIPLEWENDLANDLENS_H

#include "graleconfig.h"
#include "gravitationallens.h"
#include <vector>

namespace grale
{

class GRALE_IMPORTEXPORT WendlandLensInfo
{
public: 
	WendlandLensInfo()								{ m_heightFactor = 0; m_angularScale = 1; }
	WendlandLensInfo(double heightFactor, double angularScale, 
	                 Vector2D<double> angularPosition)				{ m_heightFactor = heightFactor; m_angularScale = angularScale; m_angularPosition = angularPosition; }

	double getHeightFactor() const							{ return m_heightFactor; }
	double getAngularScale() const							{ return m_angularScale; }
	Vector2D<double> getAngularPosition() const					{ return m_angularPosition; }
private:
	double m_heightFactor;
	double m_angularScale;
	Vector2D<double> m_angularPosition;
};

class GRALE_IMPORTEXPORT MultipleWendlandLensParams : public GravitationalLensParams
{
public:
	MultipleWendlandLensParams()							{ }
	MultipleWendlandLensParams(const std::vector<WendlandLensInfo> &phiXInfo, 
	                           const std::vector<WendlandLensInfo> &phiYInfo)		{ m_phiXInfo = phiXInfo; m_phiYInfo = phiYInfo; }

	const std::vector<WendlandLensInfo> &getPhiXInfo() const				{ return m_phiXInfo; }
	const std::vector<WendlandLensInfo> &getPhiYInfo() const				{ return m_phiYInfo; }

	void addPhiXInfo(const WendlandLensInfo &info)					{ m_phiXInfo.push_back(info); }
	void addPhiYInfo(const WendlandLensInfo &info)					{ m_phiYInfo.push_back(info); }

	bool matchDeflections(const std::vector<Vector2D<double> > &deflectionPoints,
			      const std::vector<Vector2D<double> > &deflectionAngles,
			      double angularScale);

	std::unique_ptr<GravitationalLensParams> createCopy() const;
	bool write(serut::SerializationInterface &si) const;
	bool read(serut::SerializationInterface &si);
private:
	std::vector<WendlandLensInfo> m_phiXInfo, m_phiYInfo;
};

class GRALE_IMPORTEXPORT MultipleWendlandLens : public GravitationalLens
{
public:
	MultipleWendlandLens();
	~MultipleWendlandLens();

	bool getAlphaVector(Vector2D<double> theta, Vector2D<double> *pAlpha) const;
	double getSurfaceMassDensity(Vector2D<double> theta) const ;
	bool getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *pPotentialValue) const;
	bool getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const;

	static double phiX(double r, double r2, double x, double y);
	static double phiY(double r, double r2, double x, double y);
	static double phiXX(double r, double r2, double x, double y);
	static double phiXY(double r, double r2, double x, double y);
	static double phiYY(double r, double r2, double x, double y);
	static double phiXXX(double r, double r2, double x, double y);
	static double phiXXY(double r, double r2, double x, double y);
	static double phiXYY(double r, double r2, double x, double y);
	static double phiYYY(double r, double r2, double x, double y);
protected:
	bool processParameters(const GravitationalLensParams *pLensParams);
private:
	std::vector<WendlandLensInfo> m_phiXInfo, m_phiYInfo;
};

} // end namespace

#endif // GRALE_MULTIPLEWENDLANDLENS_H
