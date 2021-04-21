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

#include "graleconfig.h"
#include "deflectionmatrix.h"
#include "gravitationallens.h"
#include <string>

// TODO
#include <iostream>

namespace grale
{

bool operator<(Vector2D<double> p, Vector2D<double> q)
{
	if (p.getX() < q.getX())
		return true;
	if (p.getX() > q.getX())
		return false;
	if (p.getY() < q.getY())
		return true;
	return false;
}

DeflectionMatrix::DeflectionMatrix()
{
	m_initialized = false;
	m_initializing = false;
}

DeflectionMatrix::~DeflectionMatrix()
{
	clear();
}

bool DeflectionMatrix::startInit()
{
	if (m_initialized)
	{
		setErrorString("Already initialized");
		return false;
	}

	if (m_initializing)
	{
		setErrorString("Alreading started initialization");
		return false;
	}

	clear(); // just to be safe

	m_initializing = true;
	return true;
}

int DeflectionMatrix::addDeflectionPoint(Vector2D<double> point)
{
	if (m_initialized)
	{
		setErrorString("Initialization procedure has already been completed");
		return -1;
	}

	if (!m_initializing)
	{
		setErrorString("Initialization procedure hasn't been started yet");
		return -1;
	}

	std::map<Vector2D<double>, int>::const_iterator it;

	if ((it = m_deflectionPointSet.find(point)) == m_deflectionPointSet.end())
	{
		int newIndex = m_deflectionPoints.size();

		m_deflectionPointSet[point] = newIndex;
		m_deflectionPoints.push_back(point);
		return newIndex;
	}

	int existingIndex = (*it).second;

	return existingIndex;
}

int DeflectionMatrix::addDerivativePoint(Vector2D<double> point)
{
	if (m_initialized)
	{
		setErrorString("Initialization procedure has already been completed");
		return -1;
	}

	if (!m_initializing)
	{
		setErrorString("Initialization procedure hasn't been started yet");
		return -1;
	}

	std::map<Vector2D<double>, int>::const_iterator it;

	if ((it = m_derivativePointSet.find(point)) == m_derivativePointSet.end())
	{
		int newIndex = m_derivativePoints.size();

		m_derivativePointSet[point] = newIndex;
		m_derivativePoints.push_back(point);
		return newIndex;
	}

	int existingIndex = (*it).second;

	return existingIndex;
}

int DeflectionMatrix::addPotentialPoint(Vector2D<double> point)
{
	if (m_initialized)
	{
		setErrorString("Initialization procedure has already been completed");
		return -1;
	}

	if (!m_initializing)
	{
		setErrorString("Initialization procedure hasn't been started yet");
		return -1;
	}

	std::map<Vector2D<double>, int>::const_iterator it;

	if ((it = m_potentialPointSet.find(point)) == m_potentialPointSet.end())
	{
		int newIndex = m_potentialPoints.size();

		m_potentialPointSet[point] = newIndex;
		m_potentialPoints.push_back(point);
		return newIndex;
	}

	int existingIndex = (*it).second;

	return existingIndex;
}

bool DeflectionMatrix::endInit(const std::vector<std::pair<std::shared_ptr<GravitationalLens>, Vector2D<double> > > &basisLenses)
{
	if (m_initialized)
	{
		setErrorString("Initialization procedure has already been completed");
		return false;
	}

	if (!m_initializing)
	{
		setErrorString("Initialization procedure hasn't been started yet");
		return false;
	}

	calculateAngularScale();

	m_numBasisFunctions = basisLenses.size();

	if (m_numBasisFunctions < 1)
	{
		clear();
		setErrorString("At least one basis function is necessary");
		return false;
	}

	// Build deflection matrix
	
	m_deflectionMatrix[0].resize(m_deflectionPoints.size());
	m_deflectionMatrix[1].resize(m_deflectionPoints.size());
	m_deflectionAngles.resize(m_deflectionPoints.size());

	for (int i = 0 ; i < m_deflectionPoints.size() ; i++)
	{
		Vector2D<double> point = m_deflectionPoints[i];

		m_deflectionMatrix[0][i].resize(m_numBasisFunctions);
		m_deflectionMatrix[1][i].resize(m_numBasisFunctions);

		for (int j = 0 ; j < m_numBasisFunctions ; j++)
		{
			Vector2D<double> deflectionAngle;
			Vector2D<double> relativePoint = point - basisLenses[j].second;
			
			if (!basisLenses[j].first->getAlphaVector(relativePoint, &deflectionAngle))
			{
				setErrorString(std::string("Couldn't calculate deflection angle: ") + basisLenses[j].first->getErrorString());
				clear();
				return false;
			}

			deflectionAngle /= m_angularScale;

			m_deflectionMatrix[0][i][j] = (float)deflectionAngle.getX();
			m_deflectionMatrix[1][i][j] = (float)deflectionAngle.getY();
		}
	}

	// Build matrix containing derivatives of the deflection angle
	
	m_derivativeMatrices[0].resize(m_derivativePoints.size());
	m_derivativeMatrices[1].resize(m_derivativePoints.size());
	m_derivativeMatrices[2].resize(m_derivativePoints.size());
	m_derivatives[0].resize(m_derivativePoints.size());
	m_derivatives[1].resize(m_derivativePoints.size());
	m_derivatives[2].resize(m_derivativePoints.size());

	for (int i = 0 ; i < m_derivativePoints.size() ; i++)
	{
		Vector2D<double> point = m_derivativePoints[i];

		m_derivativeMatrices[0][i].resize(m_numBasisFunctions);
		m_derivativeMatrices[1][i].resize(m_numBasisFunctions);
		m_derivativeMatrices[2][i].resize(m_numBasisFunctions);

		for (int j = 0 ; j < m_numBasisFunctions ; j++)
		{
			Vector2D<double> relativePoint = point - basisLenses[j].second;
			double axx, ayy, axy;

			if (!basisLenses[j].first->getAlphaVectorDerivatives(relativePoint, axx, ayy, axy))
			{
				setErrorString(std::string("Couldn't calculate deflection angle derivatives: ") + basisLenses[j].first->getErrorString());
				clear();
				return false;
			}

			m_derivativeMatrices[0][i][j] = (float)axx;
			m_derivativeMatrices[1][i][j] = (float)ayy;
			m_derivativeMatrices[2][i][j] = (float)axy;
		}
	}

	// Build potential matrix
	
	m_potentialMatrix.resize(m_potentialPoints.size());
	m_potentialValues.resize(m_potentialPoints.size());

	for (int i = 0 ; i < m_potentialPoints.size() ; i++)
	{
		Vector2D<double> point = m_potentialPoints[i];

		m_potentialMatrix[i].resize(m_numBasisFunctions);

		for (int j = 0 ; j < m_numBasisFunctions ; j++)
		{
			Vector2D<double> relativePoint = point - basisLenses[j].second;
			double potential;

			if (!basisLenses[j].first->getProjectedPotential(1.0, 1.0, relativePoint, &potential))
			{
				setErrorString(std::string("Couldn't calculate lens potential: ") + basisLenses[j].first->getErrorString());
				clear();
				return false;
			}

			potential /= (m_angularScale*m_angularScale);

			m_potentialMatrix[i][j] = (float)potential;
		}
	}

	m_initialized = true;
	m_initializing = false;

	int numBytesInMatrices = ((m_deflectionMatrix[0].size()*2 + m_derivativeMatrices[0].size()*3 + m_potentialMatrix.size())*m_numBasisFunctions*sizeof(float));

	double MB = ((double)numBytesInMatrices)/(1024.0*1024.0);

	log("Matrices need " + std::to_string(MB) + " MB");
	log("  Number of deflection points: " + std::to_string(m_deflectionMatrix[0].size()));
	log("  Number of derivative points: " + std::to_string(m_derivativeMatrices[0].size()));
	log("  Number of potential points: " + std::to_string(m_potentialMatrix.size()));
	log("  Number of basis functions: " + std::to_string(m_numBasisFunctions));

	return true;
}

void DeflectionMatrix::log(const std::string &s)
{
	std::cerr << s << std::endl;
}

void DeflectionMatrix::calculateAngularScale()
{
	double maxValue = 0;

	for (int i = 0 ; i < m_deflectionPoints.size() ; i++)
	{
		maxValue = MAX(maxValue, ABS(m_deflectionPoints[i].getX()));
		maxValue = MAX(maxValue, ABS(m_deflectionPoints[i].getY()));
	}
	for (int i = 0 ; i < m_derivativePoints.size() ; i++)
	{
		maxValue = MAX(maxValue, ABS(m_derivativePoints[i].getX()));
		maxValue = MAX(maxValue, ABS(m_derivativePoints[i].getY()));
	}
	for (int i = 0 ; i < m_potentialPoints.size() ; i++)
	{
		maxValue = MAX(maxValue, ABS(m_potentialPoints[i].getX()));
		maxValue = MAX(maxValue, ABS(m_potentialPoints[i].getY()));
	}
	
	m_angularScale = maxValue/10.0; // TODO: better criterion? Perhaps based on the actual distribution of the values?
}

bool DeflectionMatrix::calculateBasisMatrixProducts(const std::vector<float> &basisWeights, bool calcDeflection, bool calcDerivatives, bool calcPotential)
{
	if (calcDeflection)
	{
		int numPoints = m_deflectionMatrix[0].size();

		for (int i = 0 ; i < numPoints ; i++)
		{
			float x = CalculateDotProduct(&(m_deflectionMatrix[0][i][0]), &(basisWeights[0]), m_numBasisFunctions);
			float y = CalculateDotProduct(&(m_deflectionMatrix[1][i][0]), &(basisWeights[0]), m_numBasisFunctions);

			m_deflectionAngles[i] = Vector2D<float>(x, y);
		}
	}

	if (calcDerivatives)
	{
		int numPoints = m_derivativeMatrices[0].size();

		for (int i = 0 ; i < numPoints ; i++)
		{
			float axx = CalculateDotProduct(&(m_derivativeMatrices[0][i][0]), &(basisWeights[0]), m_numBasisFunctions);
			float ayy = CalculateDotProduct(&(m_derivativeMatrices[1][i][0]), &(basisWeights[0]), m_numBasisFunctions);
			float axy = CalculateDotProduct(&(m_derivativeMatrices[2][i][0]), &(basisWeights[0]), m_numBasisFunctions);

			m_derivatives[0][i] = axx;
			m_derivatives[1][i] = ayy;
			m_derivatives[2][i] = axy;
		}
	}

	if (calcPotential)
	{
		int numPoints = m_potentialMatrix.size();

		for (int i = 0 ; i < numPoints ; i++)
		{
			float potential = CalculateDotProduct(&(m_potentialMatrix[i][0]), &(basisWeights[0]), m_numBasisFunctions);

			m_potentialValues[i] = potential;
		}
	}

	return true;
}

void DeflectionMatrix::clear()
{
	m_deflectionPointSet.clear();
	m_derivativePointSet.clear();
	m_potentialPointSet.clear();

	m_deflectionPoints.clear();
	m_derivativePoints.clear();
	m_potentialPoints.clear();

	m_deflectionMatrix[0].clear();
	m_deflectionMatrix[1].clear();
	m_derivativeMatrices[0].clear();
	m_derivativeMatrices[1].clear();
	m_derivativeMatrices[2].clear();
	m_potentialMatrix.clear();

	m_deflectionAngles.clear();
	m_derivatives[0].clear();
	m_derivatives[1].clear();
	m_derivatives[2].clear();
	m_potentialValues.clear();

	m_initializing = false;
	m_initialized = false;
}

} // end namespace

