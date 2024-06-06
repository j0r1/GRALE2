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
#include "imagesbackprojector.h"
#include "imagesdataextended.h"
#include "constants.h"
#include "gravitationallens.h"
#include "lensfitnessobject.h"
#include <serut/dummyserializer.h>
#include <serut/memoryserializer.h>
#include <limits>
#include <iostream>

using namespace std;

namespace grale
{

ImagesBackProjector::ImagesBackProjector(const shared_ptr<GravitationalLens> &lens, const std::list<ImagesDataExtended *> &images,
		                         double z_d)
{
	std::list<ImagesDataExtended *>::const_iterator it;
	std::vector<ImagesDataExtended *> imagesVector(images.size());
	int i;

	for (it = images.begin(), i = 0 ; it != images.end() ; it++, i++)
		imagesVector[i] = *it;

	storeOriginalData(imagesVector);

	m_betas.resize(images.size());
	m_alphas.resize(images.size());
	m_thetas.resize(images.size());
	m_originalThetas.resize(images.size());
	m_axx.resize(images.size());
	m_ayy.resize(images.size());
	m_axy.resize(images.size());
	m_axxx.resize(images.size());
	m_ayyy.resize(images.size());
	m_axxy.resize(images.size());
	m_ayyx.resize(images.size());
	m_invMag.resize(images.size());
	m_shearComponents1.resize(images.size());
	m_convergence.resize(images.size());
	m_potential.resize(images.size());
	
	double maxx = 0, minx = 0, maxy = 0, miny = 0;
	bool extset = false;
	
	// Determine the angular scale and store some stuff

	for (it = images.begin() ; it != images.end() ; it++)
	{
		ImagesDataExtended *pImgDat = *it;
		int numImages = pImgDat->getNumberOfImages();
		
		for (int j = 0 ; j < numImages ; j++)
		{
			int numPoints = pImgDat->getNumberOfImagePoints(j);

			for (int k = 0 ; k < numPoints ; k++)
			{
				Vector2D<double> theta = pImgDat->getImagePointPosition(j, k);
				
				if (!extset)
				{
					minx = theta.getX();
					maxx = theta.getX();
					miny = theta.getY();
					maxy = theta.getY();
					extset = true;
				}
				else
				{
					minx = MIN(minx,theta.getX());
					maxx = MAX(maxx,theta.getX());
					miny = MIN(miny,theta.getY());
					maxy = MAX(maxy,theta.getY());
				}
			}			
		}
	}
	
	m_angularScale = MIN(maxx-minx, maxy-miny)/15.0; // Let's use something different than in deflectionmatrix
	
	// Store thetas

	for (int s = 0 ; s < imagesVector.size() ;  s++)
	{
		ImagesDataExtended *pImgDat = imagesVector[s];

		m_thetas[s].resize(m_numTotalPoints[s]);
		m_originalThetas[s].resize(m_numTotalPoints[s]);
		int numImages = pImgDat->getNumberOfImages();
		int pos = 0;

		for (int j = 0 ; j < numImages ; j++)
		{
			int numPoints = pImgDat->getNumberOfImagePoints(j);

			for (int k = 0 ; k < numPoints ; k++, pos++)
			{
				Vector2D<double> theta = pImgDat->getImagePointPosition(j,k);
				
				m_originalThetas[s][pos] = theta;

				theta /= m_angularScale;
				
				m_thetas[s][pos] = Vector2D<float>(theta.getX(), theta.getY());
			}			
		}
	}

	m_pLens = lens;
	m_zd = z_d;
}

ImagesBackProjector::~ImagesBackProjector()
{
}

float ImagesBackProjector::getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> scaledBeta) const
{
	double td = 0;
	Vector2D<double> theta = m_originalThetas[sourceNumber][m_offsets[sourceNumber][imageNumber] + pointNumber];
	Vector2D<double> beta(scaledBeta.getX(),scaledBeta.getY());
	
	beta *= m_angularScale;

	m_pLens->getTimeDelay(m_zd, 1.0, m_distanceFractions[sourceNumber], theta, beta, &td);
	td /= (60.0*60.0*24.0); // convert to days
	return (float)td;
}

void ImagesBackProjector::checkBetas(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_betas[sourceNumber].size() == numPoints)
		return;

	m_betas[sourceNumber].resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
	{
		Vector2D<double> beta { numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN() };

		m_pLens->traceTheta(1.0, m_distanceFractions[sourceNumber], m_originalThetas[sourceNumber][i], &beta);

		beta /= m_angularScale;
		m_betas[sourceNumber][i] = Vector2D<float>(beta.getX(), beta.getY());
	}
}

void ImagesBackProjector::checkAlphas(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_alphas[sourceNumber].size() == numPoints)
		return;

	m_alphas[sourceNumber].resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
	{
		Vector2D<double> alpha { numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN() };

		m_pLens->getAlphaVector(m_originalThetas[sourceNumber][i], &alpha);
		alpha *= m_distanceFractions[sourceNumber];
		alpha /= m_angularScale;
		m_alphas[sourceNumber][i] = Vector2D<float>(alpha.getX(), alpha.getY());
	}
}

void ImagesBackProjector::checkDerivatives(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_axx[sourceNumber].size() == numPoints)
		return;

	m_axx[sourceNumber].resize(numPoints);
	m_ayy[sourceNumber].resize(numPoints);
	m_axy[sourceNumber].resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
	{
		double axx = numeric_limits<double>::quiet_NaN();
		double ayy = numeric_limits<double>::quiet_NaN();
		double axy = numeric_limits<double>::quiet_NaN();

		m_pLens->getAlphaVectorDerivatives(m_originalThetas[sourceNumber][i], axx, ayy, axy);

		axx *= m_distanceFractions[sourceNumber];
		axy *= m_distanceFractions[sourceNumber];
		ayy *= m_distanceFractions[sourceNumber];

		m_axx[sourceNumber][i] = (float)axx;
		m_ayy[sourceNumber][i] = (float)ayy;
		m_axy[sourceNumber][i] = (float)axy;
	}
}

void ImagesBackProjector::checkSecondDerivatives(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_axxx[sourceNumber].size() == numPoints)
		return;

	m_axxx[sourceNumber].resize(numPoints);
	m_ayyy[sourceNumber].resize(numPoints);
	m_axxy[sourceNumber].resize(numPoints);
	m_ayyx[sourceNumber].resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
	{
		double axxx = numeric_limits<double>::quiet_NaN();
		double ayyy = numeric_limits<double>::quiet_NaN();
		double axxy = numeric_limits<double>::quiet_NaN();
		double ayyx = numeric_limits<double>::quiet_NaN();

		m_pLens->getAlphaVectorSecondDerivatives(m_originalThetas[sourceNumber][i], axxx, ayyy, axxy, ayyx);

		axxx *= m_distanceFractions[sourceNumber];
		ayyy *= m_distanceFractions[sourceNumber];
		axxy *= m_distanceFractions[sourceNumber];
		ayyx *= m_distanceFractions[sourceNumber];

		m_axxx[sourceNumber][i] = (float)(axxx*m_angularScale);
		m_ayyy[sourceNumber][i] = (float)(ayyy*m_angularScale);
		m_axxy[sourceNumber][i] = (float)(axxy*m_angularScale);
		m_ayyx[sourceNumber][i] = (float)(ayyx*m_angularScale);
	}
}

void ImagesBackProjector::checkInvMag(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_invMag[sourceNumber].size() == numPoints)
		return;

	checkDerivatives(sourceNumber);

	m_invMag[sourceNumber].resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
		m_invMag[sourceNumber][i] = (float)((1.0-(double)m_axx[sourceNumber][i])*(1.0-(double)m_ayy[sourceNumber][i])-m_axy[sourceNumber][i]*m_axy[sourceNumber][i]);
}

void ImagesBackProjector::checkShear(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_shearComponents1[sourceNumber].size() == numPoints)
		return;

	checkDerivatives(sourceNumber);

	m_shearComponents1[sourceNumber].resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
		m_shearComponents1[sourceNumber][i] = (float)(0.5*((double)m_axx[sourceNumber][i] - (double)m_ayy[sourceNumber][i]));
}

void ImagesBackProjector::checkConvergence(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_convergence[sourceNumber].size() == numPoints)
		return;

	checkDerivatives(sourceNumber);

	m_convergence[sourceNumber].resize(numPoints);

	for (int i = 0 ; i < numPoints ; i++)
		m_convergence[sourceNumber][i] = (float)(0.5*((double)m_axx[sourceNumber][i] + (double)m_ayy[sourceNumber][i]));
}

void ImagesBackProjector::checkPotential(int sourceNumber) const
{
	int numPoints = m_thetas[sourceNumber].size();

	if (m_potential[sourceNumber].size() == numPoints)
		return;

	m_potential[sourceNumber].resize(numPoints);
	double potentialScale = m_angularScale*m_angularScale;

	for (int i = 0 ; i < numPoints ; i++)
	{
		double potential = 0;
		m_pLens->getProjectedPotential(1.0, m_distanceFractions[sourceNumber], m_originalThetas[sourceNumber][i], &potential);
		m_potential[sourceNumber][i] = (float)(potential/potentialScale);
	}
}

double ImagesBackProjector::getLensDistance() const
{ 
	return m_pLens->getLensDistance(); 
}

} // end namespace

