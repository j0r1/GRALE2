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
#include "backprojectmatrixnew.h"
#include "gravitationallens.h"
#include "imagesdataextended.h"
#include "constants.h"
#include <errut/booltype.h>

// TODO
#include <iostream>

#include "debugnew.h"

using namespace std;
using namespace errut;

namespace grale
{

BackProjectMatrixNew::BackProjectMatrixNew()
{
	m_init = false;
	m_initializing = false;
}

BackProjectMatrixNew::~BackProjectMatrixNew()
{
}

bool BackProjectMatrixNew::startInit(double z_d, double D_d, DeflectionMatrix *pDeflectionMatrix, 
		                     const std::vector<ImagesDataExtended *> &images,
		                     const std::vector<bool> &useDeflections, 
				     const std::vector<bool> &useDerivatives, 
				     const std::vector<bool> &usePotentials, 
				     const GravitationalLens *pBaseLens,
		       		     bool storeOriginalIntensities,
		       		     bool storeOriginalTimeDelays,
		       		     bool storeOriginalShearInfo,
					 const GravitationalLens *pSheetLens)
{
	if (m_init)
	{
		setErrorString("Already inialized");
		return false;
	}
	if (m_initializing)
	{
		setErrorString("Already started initialization");
		return false;
	}

	m_initializing = true;
	m_pDeflectionMatrix = pDeflectionMatrix;

	storeOriginalData(images, storeOriginalIntensities, storeOriginalTimeDelays, storeOriginalShearInfo);

	m_deflectionIndices.resize(images.size());
	m_derivativeIndices.resize(images.size());
	m_potentialIndices.resize(images.size());

	m_originalPoints.resize(images.size());
	m_thetas.resize(images.size());
	m_subDeflectionAngles.resize(images.size());
	m_subDeflectionDerivatives[0].resize(images.size());
	m_subDeflectionDerivatives[1].resize(images.size());
	m_subDeflectionDerivatives[2].resize(images.size());
	m_subPotentialValues.resize(images.size());

	m_betas.resize(images.size());
	m_alphas.resize(images.size());
	m_axx.resize(images.size());
	m_ayy.resize(images.size());
	m_axy.resize(images.size());
	m_potentials.resize(images.size());
	m_inverseMagnifications.resize(images.size());
	m_shearComponent1.resize(images.size());
	m_convergence.resize(images.size());

	if (pBaseLens)
	{
		m_baseAlphas.resize(images.size());
		m_baseAxx.resize(images.size());
		m_baseAyy.resize(images.size());
		m_baseAxy.resize(images.size());
		m_basePotentials.resize(images.size());
		m_useBaseLens = true;
	}
	else
		m_useBaseLens = false;

	if (pSheetLens)
	{
		m_useMassSheet = true;
		m_sheetAlphas.resize(images.size());
		m_sheetAxx.resize(images.size());
		m_sheetAyy.resize(images.size());
		m_sheetAxy.resize(images.size());
		m_sheetPotentials.resize(images.size());
	}
	else
		m_useMassSheet = false;

	int maxPoints = 0;

	for (int s = 0 ; s < images.size() ; s++)
	{
		ImagesDataExtended *pImgDat = images[s];
		int numImages = pImgDat->getNumberOfImages();
		int offset = m_numTotalPoints[s];

		m_originalPoints[s].resize(offset); // we'll store all the points since they may be necessary if a base lens is used
		m_thetas[s].resize(offset); // These are also necessary for time delay calculations

		if (useDeflections[s])
		{
			m_deflectionIndices[s].resize(offset);

			for (int i = 0 ; i < m_deflectionIndices[s].size() ; i++)
				m_deflectionIndices[s][i] = -1;
			m_subDeflectionAngles[s].resize(offset);
			m_betas[s].resize(offset);
			m_alphas[s].resize(offset);
			if (m_useBaseLens)
				m_baseAlphas[s].resize(offset);
			if (m_useMassSheet)
				m_sheetAlphas[s].resize(offset);
		}

		if (useDerivatives[s])
		{
			m_derivativeIndices[s].resize(offset);
			for (int i = 0 ; i < m_derivativeIndices[s].size() ; i++)
				m_derivativeIndices[s][i] = -1;
			m_subDeflectionDerivatives[0][s].resize(offset);
			m_subDeflectionDerivatives[1][s].resize(offset);
			m_subDeflectionDerivatives[2][s].resize(offset);
			m_axx[s].resize(offset);
			m_ayy[s].resize(offset);
			m_axy[s].resize(offset);
			if (m_useBaseLens)
			{
				m_baseAxx[s].resize(offset);
				m_baseAyy[s].resize(offset);
				m_baseAxy[s].resize(offset);
			}
			if (m_useMassSheet)
			{
				m_sheetAxx[s].resize(offset);
				m_sheetAyy[s].resize(offset);
				m_sheetAxy[s].resize(offset);
			}

			m_inverseMagnifications[s].resize(offset);
			m_shearComponent1[s].resize(offset);
			m_convergence[s].resize(offset);
		}

		if (pImgDat->hasTimeDelays() && usePotentials[s])
		{
			m_potentialIndices[s].resize(offset);
			for (int i = 0 ; i < m_potentialIndices[s].size() ; i++)
				m_potentialIndices[s][i] = -1;
			m_subPotentialValues[s].resize(offset);
			m_potentials[s].resize(offset);
			if (m_useBaseLens)
				m_basePotentials[s].resize(offset);
			if (m_useMassSheet)
				m_sheetPotentials[s].resize(offset);
		}

		if (offset > maxPoints)
			maxPoints = offset;
	}

	for (int s = 0 ; s < images.size() ; s++)
	{
		ImagesDataExtended *pImgDat = images[s];
		int numImages = pImgDat->getNumberOfImages();
		int offset = 0;

		for (int i = 0 ; i < numImages ; i++)
		{
			int numPoints = pImgDat->getNumberOfImagePoints(i);

			for (int p = 0 ; p < numPoints ; p++, offset++)
			{
				Vector2D<double> point = pImgDat->getImagePointPosition(i, p);
				int pointIndex;

				m_originalPoints[s][offset] = point;

				if (useDeflections[s])
				{
					if ((pointIndex = m_pDeflectionMatrix->addDeflectionPoint(point)) < 0)
					{
						setErrorString(std::string("Couldn't add deflection point: ") + m_pDeflectionMatrix->getErrorString());
						return false;
					}

					m_deflectionIndices[s][offset] = pointIndex;
				}
				if (useDerivatives[s])
				{
					if ((pointIndex = m_pDeflectionMatrix->addDerivativePoint(point)) < 0)
					{
						setErrorString(std::string("Couldn't add derivative point: ") + m_pDeflectionMatrix->getErrorString());
						return false;
					}

					m_derivativeIndices[s][offset] = pointIndex;
				}

			}
		}

		if (usePotentials[s])
		{
			int numTimeDelays = pImgDat->getNumberOfTimeDelays();

			for (int i = 0 ; i < numTimeDelays ; i++)
			{
				int imageNumber, pointNumber;
				double delay;

				pImgDat->getTimeDelay(i, &imageNumber, &pointNumber, &delay);

				int pointIndex;
				int offset = m_offsets[s][imageNumber] + pointNumber;
				Vector2D<double> point = m_originalPoints[s][offset];
				
				if ((pointIndex = m_pDeflectionMatrix->addPotentialPoint(point)) < 0)
				{
					setErrorString(std::string("Couldn't add potential point: ") + m_pDeflectionMatrix->getErrorString());
					return false;
				}
				
				m_potentialIndices[s][offset] = pointIndex;
			}
		}
	}

	auto &originalPoints = m_originalPoints;
	auto &distanceFractions = m_distanceFractions;
	auto getLensProperties = [&originalPoints, &distanceFractions](const GravitationalLens *pLens,
								vector<vector<Vector2Dd>> &alphasUnscaled, 
								vector<vector<double>> &potentialsUnscaled,
								vector<vector<Vector2Df>> &alphas,
								vector<vector<float>> &potentials,
								vector<vector<float>> &axxs,
								vector<vector<float>> &ayys,
								vector<vector<float>> &axys) -> bool_t
	{
		alphasUnscaled.resize(originalPoints.size());
		potentialsUnscaled.resize(originalPoints.size());

		for (int s = 0 ; s < originalPoints.size() ; s++)
		{
			int numPoints;

			if ((numPoints = alphas[s].size()) > 0)
			{
				alphasUnscaled[s].resize(numPoints);

				for (int i = 0 ; i < numPoints ; i++)
				{
					Vector2D<double> point = originalPoints[s][i];
					Vector2D<double> alpha;

					if (!pLens->getAlphaVector(point, &alpha))
						return "Unable to calculate base lens bending angle: " + pLens->getErrorString();

					alphasUnscaled[s][i] = ((double)distanceFractions[s]) * alpha;
				}
			}
			if ((numPoints = axxs[s].size()) > 0)
			{
				for (int i = 0 ; i < numPoints ; i++)
				{
					Vector2D<double> point = originalPoints[s][i];
					double axx, ayy, axy;

					if (!pLens->getAlphaVectorDerivatives(point, axx, ayy, axy))
						return "Unable to calculate base lens deflection angle derivatives: " + pLens->getErrorString();

					axx *= (double)distanceFractions[s];
					ayy *= (double)distanceFractions[s];
					axy *= (double)distanceFractions[s];

					axxs[s][i] = (float)axx;
					ayys[s][i] = (float)ayy;
					axys[s][i] = (float)axy;
				}
			}
			if ((numPoints = potentials[s].size()) > 0)
			{
				potentialsUnscaled[s].resize(numPoints);

				for (int i = 0 ; i < numPoints ; i++)
				{
					Vector2D<double> point = originalPoints[s][i];
					double potential;

					if (!pLens->getProjectedPotential(1.0, distanceFractions[s], point, &potential))
						return "Unable to calculate base lens potential: " + pLens->getErrorString();
					
					potentialsUnscaled[s][i] = potential;
				}
			}
		}
		return true;
	};

	if (m_useBaseLens)
	{
		bool_t r = getLensProperties(pBaseLens, m_baseAlphasUnscaled, m_basePotentialsUnscaled, 
				                     m_baseAlphas, m_basePotentials, m_baseAxx, m_baseAyy, m_baseAxy);
		if (!r)
		{
			setErrorString("Can't get lens properties for base lens: " + r.getErrorString());
			return false;
		}
	}
	if (m_useMassSheet)
	{
		bool_t r = getLensProperties(pSheetLens, m_sheetAlphasUnscaled, m_sheetPotentialsUnscaled,
				                     m_sheetAlphas, m_sheetPotentials, m_sheetAxx, m_sheetAyy, m_sheetAxy);

		if (!r)
		{
			setErrorString("Can't get lens properties for sheet lens: " + r.getErrorString());
			return false;
		}
	}

	m_tmpBuffer.resize(maxPoints*2); // *2 since it can be about Vector2D<float>
	m_oneVector.resize(maxPoints);

	for (int i = 0 ; i < maxPoints ; i++)
		m_oneVector[i] = 1;

	m_trueFlags.resize(images.size());
	for (int i = 0 ; i < images.size() ; i++)
		m_trueFlags[i] = true;

	m_Dd = D_d;
	m_zd = z_d;

	return true;
}

bool BackProjectMatrixNew::endInit() // m_pDeflectionMatrix->endInit() must be called before this is called (to know the angular scale)
{
	if (m_init)
	{
		setErrorString("Already initialized");
		return false;
	}
	if (!m_initializing)
	{
		setErrorString("Initialization hasn't been started yet");
		return false;
	}

	double angularScale = m_pDeflectionMatrix->getAngularScale();

	for (int s = 0 ; s < m_originalPoints.size() ; s++)
	{
		for (int i = 0 ; i < m_originalPoints[s].size() ; i++)
		{
			Vector2D<double> point = m_originalPoints[s][i];

			m_thetas[s][i] = Vector2D<float>((float)(point.getX()/angularScale), (float)(point.getY()/angularScale));
		}

		if (m_useBaseLens)
		{
			int numPoints;

			if ((numPoints = m_baseAlphas[s].size()) > 0)
			{
				for (int i = 0 ; i < numPoints ; i++)
					m_baseAlphas[s][i] = Vector2D<float>(m_baseAlphasUnscaled[s][i].getX()/angularScale, m_baseAlphasUnscaled[s][i].getY()/angularScale);
			}

			if ((numPoints = m_basePotentials[s].size()) > 0)
			{
				for (int i = 0 ; i < numPoints ; i++)
					m_basePotentials[s][i] = m_basePotentialsUnscaled[s][i]/(angularScale*angularScale);
			}
		}

		if (m_useMassSheet)
		{
			int numPoints;

			if ((numPoints = m_sheetAlphas[s].size()) > 0)
			{
				for (int i = 0 ; i < numPoints ; i++)
					m_sheetAlphas[s][i] = Vector2D<float>(m_sheetAlphasUnscaled[s][i].getX()/angularScale, m_sheetAlphasUnscaled[s][i].getY()/angularScale);
			}
			if ((numPoints = m_sheetPotentials[s].size()) > 0)
			{
				for (int i = 0 ; i < numPoints ; i++)
					m_sheetPotentials[s][i] = m_sheetPotentialsUnscaled[s][i]/(angularScale*angularScale);
			}
		}
	}

	m_originalPoints.clear();
	m_baseAlphasUnscaled.clear();
	m_basePotentialsUnscaled.clear();

	m_timeDelayScale = (float)((1.0+m_zd)*m_Dd*angularScale*angularScale/(SPEED_C*60.0*60.0*24.0));
	m_massSheetScale = SPEED_C*SPEED_C/(4.0*CONST_PI*CONST_G*m_Dd);
	m_angularScale = angularScale;

	m_initializing = false;
	m_init = true;
	return true;
}

void BackProjectMatrixNew::storeDeflectionMatrixResults()
{
	int numSources = m_offsets.size();

	for (int s = 0 ; s < numSources ; s++)
	{
		float distanceFraction = m_distanceFractions[s];
		int numPoints;

		numPoints = m_deflectionIndices[s].size();

		for (int i = 0 ; i < numPoints ; i++)
		{
			int index = m_deflectionIndices[s][i];

			m_subDeflectionAngles[s][i] = m_pDeflectionMatrix->getDeflectionAngle(index) * distanceFraction;
		}

		numPoints = m_derivativeIndices[s].size();

		for (int i = 0 ; i < numPoints ; i++)
		{
			int index = m_derivativeIndices[s][i];
			float axx, ayy, axy;

			m_pDeflectionMatrix->getDeflectionDerivatives(index, &axx, &ayy, &axy);
			m_subDeflectionDerivatives[0][s][i] = axx * distanceFraction;
			m_subDeflectionDerivatives[1][s][i] = ayy * distanceFraction;
			m_subDeflectionDerivatives[2][s][i] = axy * distanceFraction;
		}

		numPoints = m_potentialIndices[s].size();

		for (int i = 0 ; i < numPoints ; i++)
		{
			int index = m_potentialIndices[s][i];

			if (index >= 0)
				m_subPotentialValues[s][i] = m_pDeflectionMatrix->getPotential(index) * distanceFraction;
			else
				m_subPotentialValues[s][i] = 0;
		}
	}
}

void BackProjectMatrixNew::calculate(float scaleFactor, float massSheetFactor)
{
	int numSources = m_offsets.size();

	for (int s = 0 ; s < numSources ; s++)
	{
		int numPoints;

		if ((numPoints = m_deflectionIndices[s].size()) > 0)
		{
			MultiplyVector((float *)&(m_alphas[s][0]), (float *)&(m_subDeflectionAngles[s][0]), scaleFactor, numPoints*2);

			if (m_useBaseLens)
				AddVector((float *)&(m_alphas[s][0]), (float *)&(m_baseAlphas[s][0]), numPoints*2);
			if (massSheetFactor != 0 && m_useMassSheet)
				CalculateAddProductC((float *)&(m_alphas[s][0]), (float *)&(m_sheetAlphas[s][0]), massSheetFactor, numPoints*2, &(m_tmpBuffer[0]));

			SubVector((float *)&(m_betas[s][0]), (float *)&(m_thetas[s][0]), (float *)&(m_alphas[s][0]), numPoints*2);
		}

		if ((numPoints = m_derivativeIndices[s].size()) > 0)
		{
			MultiplyVector(&(m_axx[s][0]), &(m_subDeflectionDerivatives[0][s][0]), scaleFactor, numPoints);
			MultiplyVector(&(m_ayy[s][0]), &(m_subDeflectionDerivatives[1][s][0]), scaleFactor, numPoints);
			MultiplyVector(&(m_axy[s][0]), &(m_subDeflectionDerivatives[2][s][0]), scaleFactor, numPoints);
			if (m_useBaseLens)
			{
				AddVector(&(m_axx[s][0]), &(m_baseAxx[s][0]), numPoints);
				AddVector(&(m_ayy[s][0]), &(m_baseAyy[s][0]), numPoints);
				AddVector(&(m_axy[s][0]), &(m_baseAxy[s][0]), numPoints);
			}
			if (massSheetFactor != 0 && m_useMassSheet)
			{
				//AddVector(&(m_axx[s][0]), massSheetFactor*m_distanceFractions[s], numPoints);
				//AddVector(&(m_ayy[s][0]), massSheetFactor*m_distanceFractions[s], numPoints);
				// axy is zero for mass sheet
				CalculateAddProductC(&(m_axx[s][0]), (float *)&(m_sheetAxx[s][0]), massSheetFactor, numPoints, &(m_tmpBuffer[0]));
				CalculateAddProductC(&(m_ayy[s][0]), (float *)&(m_sheetAyy[s][0]), massSheetFactor, numPoints, &(m_tmpBuffer[0]));
				CalculateAddProductC(&(m_axy[s][0]), (float *)&(m_sheetAxy[s][0]), massSheetFactor, numPoints, &(m_tmpBuffer[0]));
			}
		}

		if ((numPoints = m_potentialIndices[s].size()) > 0)
		{
			MultiplyVector(&(m_potentials[s][0]), &(m_subPotentialValues[s][0]), scaleFactor, numPoints);
			if (m_useBaseLens)
				AddVector(&(m_potentials[s][0]), &(m_basePotentials[s][0]), numPoints);

			if (massSheetFactor != 0 && m_useMassSheet)
				CalculateAddProductC(&(m_potentials[s][0]), &(m_sheetPotentials[s][0]), massSheetFactor, numPoints, &(m_tmpBuffer[0]));
		}
	}
}

void BackProjectMatrixNew::calculateInverseMagnifications(const std::vector<bool> &sourceMask)
{
	int numSources = m_offsets.size();

	for (int s = 0 ; s < numSources ; s++)
	{
		int numPoints = m_derivativeIndices[s].size();

		if (sourceMask[s] && (numPoints > 0))
		{
			SubVector(&(m_inverseMagnifications[s][0]), &(m_oneVector[0]), &(m_axx[s][0]), numPoints);
			SubVector(&(m_inverseMagnifications[s][0]), &(m_ayy[s][0]), numPoints);
			MultiplyVector(&(m_tmpBuffer[0]), &(m_axx[s][0]), &(m_ayy[s][0]), numPoints);
			AddVector(&(m_inverseMagnifications[s][0]), &(m_tmpBuffer[0]), numPoints);
			SquareVector(&(m_tmpBuffer[0]), &(m_axy[s][0]), numPoints);
			SubVector(&(m_inverseMagnifications[s][0]), &(m_tmpBuffer[0]), numPoints);
		}
	}
}

void BackProjectMatrixNew::calculateShearComponents(const std::vector<bool> &sourceMask)
{
	int numSources = m_offsets.size();

	for (int s = 0 ; s < numSources ; s++)
	{
		int numPoints = m_derivativeIndices[s].size();

		if (sourceMask[s] && (numPoints > 0))
		{
			SubVector(&(m_shearComponent1[s][0]), &(m_axx[s][0]), &(m_ayy[s][0]), numPoints);
			MultiplyVector(&(m_shearComponent1[s][0]), 0.5f, numPoints);
		}
	}
}

void BackProjectMatrixNew::calculateConvergence(const std::vector<bool> &sourceMask)
{
	int numSources = m_offsets.size();

	for (int s = 0 ; s < numSources ; s++)
	{
		int numPoints = m_derivativeIndices[s].size();

		if (sourceMask[s] && (numPoints > 0))
		{
			AddVector(&(m_convergence[s][0]), &(m_axx[s][0]), &(m_ayy[s][0]), numPoints);
			MultiplyVector(&(m_convergence[s][0]), 0.5f, numPoints);
		}
	}
}

} // end namespace

