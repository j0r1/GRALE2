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

#define NOMINMAX
#include "galensmodule.h" // includes windows.h, causes problems with std::min and std::max

#include <iostream>

#include "debugnew.h"

namespace grale
{

ImagesBackProjector::ImagesBackProjector(GravitationalLens &lens, const std::list<ImagesDataExtended *> &images, 
		                         double z_d, bool copyLens)
{
	std::list<ImagesDataExtended *>::const_iterator it;
	std::vector<ImagesDataExtended *> imagesVector(images.size());
	int i;

	for (it = images.begin(), i = 0 ; it != images.end() ; it++, i++)
		imagesVector[i] = *it;

	storeOriginalData(imagesVector, true, true, true);

	m_betas.resize(images.size());
	m_alphas.resize(images.size());
	m_thetas.resize(images.size());
	m_originalThetas.resize(images.size());
	m_axx.resize(images.size());
	m_ayy.resize(images.size());
	m_axy.resize(images.size());
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

	if (copyLens)
	{
		m_pLens = lens.createCopy();
		m_deleteLens = true;
	}
	else
	{
		m_pLens = &lens;
		m_deleteLens = false;
	}
	m_zd = z_d;
}

ImagesBackProjector::~ImagesBackProjector()
{
	if (m_deleteLens)
		delete m_pLens;
}

bool ImagesBackProjector::write(const std::string &fname, double angularunit, bool magnifyFlux) const
{
	FILE *f = fopen(fname.c_str(), "wt");
	if (f == 0)
	{
		setErrorString("Can't open specified file");
		return false;
	}

	for (int s = 0 ; s < m_numPoints.size() ; s++)
	{
		for (int i = 0 ; i < m_numPoints[s].size() ; i++)
		{
			int numPoints = m_numPoints[s][i];

			for (int p = 0 ; p < numPoints ; p++)
			{
				Vector2D<float> beta = getBetas(s, i)[p];
				beta *= (m_angularScale/angularunit);

				double height = 0;

				if (hasOriginalIntensities(s))
				{
					double intens = getOriginalIntensities(s, i)[p];

					if (magnifyFlux)
						height = intens*ABS(getInverseMagnifications(s, i)[p])*getIntensityScale();
					else
						height = intens*getIntensityScale();
				}

				fprintf(f, "%g %g %g\n", (double)beta.getX(), (double)beta.getY(), (double)height);
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\n\n");
	}

	fclose(f);
	return true;
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
		Vector2D<double> beta;

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
		Vector2D<double> alpha(0,0);

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
		double axx, ayy, axy;

		m_pLens->getAlphaVectorDerivatives(m_originalThetas[sourceNumber][i], axx, ayy, axy);

		axx *= m_distanceFractions[sourceNumber];
		axy *= m_distanceFractions[sourceNumber];
		ayy *= m_distanceFractions[sourceNumber];

		m_axx[sourceNumber][i] = (float)axx;
		m_ayy[sourceNumber][i] = (float)ayy;
		m_axy[sourceNumber][i] = (float)axy;
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

} // end namespace

