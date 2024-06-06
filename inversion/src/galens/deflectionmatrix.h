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

#ifndef GRALE_DEFLECTIONMATRIX_H

#define GRALE_DEFLECTIONMATRIX_H

#include "graleconfig.h"
#include "vector2d.h"
#include <errut/errorbase.h>
#include <vector>
#include <map>
#include <memory>
#include <cassert>

namespace grale
{

class GravitationalLens;

class GRALE_IMPORTEXPORT DeflectionMatrix : public errut::ErrorBase
{
public:
	DeflectionMatrix();
	~DeflectionMatrix();

	bool startInit();
	int addDeflectionPoint(Vector2D<double> point);
	int addDerivativePoint(Vector2D<double> point);
	int addPotentialPoint(Vector2D<double> point);
	int addSecondDerivativePoint(Vector2D<double> point);
	bool endInit(const std::vector<std::pair<std::shared_ptr<GravitationalLens>, Vector2D<double> > > &basisLenses);

	bool calculateBasisMatrixProducts(const std::vector<float> &basisWeights, bool calcDeflection, bool calcDerivatives, bool calcPotential,
	                                  bool calcSecondDerivatives);

	double getAngularScale() const									{ return m_angularScale; }

	Vector2D<float> getDeflectionAngle(int index) const;
	void getDeflectionDerivatives(int index, float *pAxx, float *pAyy, float *pAxy) const;
	float getPotential(int index) const;
	void getSecondDeflectionDerivatives(int index, float *pAxxx, float *pAyyy, float *pAxxy, float *pAyyx) const;
protected:
	virtual void log(const std::string &s);
private:
	void calculateAngularScale();
	void clear();

	bool m_initializing;
	bool m_initialized;

	// Used during initialization, to determine for which points which values need to be
	// calculated.

	std::map<Vector2D<double>, int> m_deflectionPointSet;
	std::map<Vector2D<double>, int> m_derivativePointSet;
	std::map<Vector2D<double>, int> m_potentialPointSet;
	std::map<Vector2D<double>, int> m_secondDerivativePointSet;

	// The actual points for which the deflections are being calculated, we'll keep these
	// for debugging purposes

	std::vector<Vector2D<double> > m_deflectionPoints;
	std::vector<Vector2D<double> > m_derivativePoints;
	std::vector<Vector2D<double> > m_potentialPoints;
	std::vector<Vector2D<double> > m_secondDerivativePoints;

	double m_angularScale;
	int m_numBasisFunctions;

	std::vector<std::vector<float > > m_deflectionMatrix[2]; // stores for each point in the lens plane the deflection angle by each basis function
	std::vector<std::vector<float> > m_derivativeMatrices[3]; // same, but for the derivatives of the deflection angle
	std::vector<std::vector<float> > m_potentialMatrix; // same, but for the lens potential
	std::vector<std::vector<float> > m_secondDerivativeMatrices[4]; // same, for second derivs of alpha

	std::vector<Vector2D<float> > m_deflectionAngles; // resulting deflection angles, for a specific set of weights
	std::vector<float> m_derivatives[3]; // resulting derivatives
	std::vector<float> m_potentialValues; // resulting values of the lens potential
	std::vector<float> m_secondDerivatives[4]; // resulting second derivatives
};

inline Vector2D<float> DeflectionMatrix::getDeflectionAngle(int index) const
{
	assert(index >= 0 && (size_t)index < m_deflectionAngles.size());
	return m_deflectionAngles[index];
}

inline void DeflectionMatrix::getDeflectionDerivatives(int index, float *pAxx, float *pAyy, float *pAxy) const
{
	assert(index >= 0 && (size_t)index < m_derivatives[0].size() &&
						 (size_t)index < m_derivatives[1].size() &&
						 (size_t)index < m_derivatives[2].size());
	assert(pAxx && pAyy && pAxy);

	*pAxx = m_derivatives[0][index];
	*pAyy = m_derivatives[1][index];
	*pAxy = m_derivatives[2][index];
}

inline float DeflectionMatrix::getPotential(int index) const
{
	assert(index >= 0 && (size_t)index < m_potentialValues.size());
	return m_potentialValues[index];
}

inline void DeflectionMatrix::getSecondDeflectionDerivatives(int index, float *pAxxx, float *pAyyy, float *pAxxy, float *pAyyx) const
{
	assert(index >= 0 && (size_t)index < m_secondDerivatives[0].size() &&
						 (size_t)index < m_secondDerivatives[1].size() &&
						 (size_t)index < m_secondDerivatives[2].size() &&
						 (size_t)index < m_secondDerivatives[3].size());
	assert(pAxxx && pAyyy && pAxxy && pAyyx);

	*pAxxx = m_secondDerivatives[0][index];
	*pAyyy = m_secondDerivatives[1][index];
	*pAxxy = m_secondDerivatives[2][index];
	*pAyyx = m_secondDerivatives[3][index];
}

} // end namespace

#endif // GRALE_DEFLECTIONMATRIX_H

