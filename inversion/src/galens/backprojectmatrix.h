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

#pragma once

/**
 * \file backprojectmatrix.h
 */

#include "graleconfig.h"
#include "deflectionmatrix.h"
#include "projectedimagesinterface.h"

namespace grale
{

class ImagesDataExtended;
class GravitationalLens;

/** Implements the ProjectedImagesInterface interface, is only meant to be used inside
 *  the genetic algorithm for lens inversion. */
class GRALE_IMPORTEXPORT BackProjectMatrix : public ProjectedImagesInterface, public errut::ErrorBase
{
public:
	BackProjectMatrix();
	~BackProjectMatrix();

	bool startInit(double z_d, double D_d, DeflectionMatrix *pDeflectionMatrix, 
	               const std::vector<ImagesDataExtended *> &images,
	               const std::vector<bool> &useDeflections, 
		       const std::vector<bool> &useDerivatives,
		       const std::vector<bool> &usePotentials,
			   const std::vector<bool> &useSecondDerivs,
		       const GravitationalLens *pBaseLens,
		       const GravitationalLens *pSheetLens);
	bool endInit();

	void storeDeflectionMatrixResults();
	void calculate(float scaleFactor, float massSheetFactor = 0);
	void calculateInverseMagnifications()									{ calculateInverseMagnifications(m_trueFlags); }
	void calculateShearComponents()										{ calculateShearComponents(m_trueFlags); }
	void calculateConvergence()										{ calculateConvergence(m_trueFlags); }
	void calculateInverseMagnifications(const std::vector<bool> &sourceMask);
	void calculateShearComponents(const std::vector<bool> &sourceMask);
	void calculateConvergence(const std::vector<bool> &sourceMask);

	double getAngularScale() const										{ return m_angularScale; }
	double getMassSheetScale() const									{ return m_massSheetScale; }

	double getLensDistance() const										{ return m_Dd; }
	double getLensRedshift() const										{ return m_zd; }
	const Vector2D<float> *getBetas(int sourceNumber) const							{ return &(m_betas[sourceNumber][0]); }
	const Vector2D<float> *getBetas(int sourceNumber, int imageNumber) const 				{ return &(m_betas[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const Vector2D<float> *getAlphas(int sourceNumber) const							{ return &(m_alphas[sourceNumber][0]); }
	const Vector2D<float> *getAlphas(int sourceNumber, int imageNumber) const 				{ return &(m_alphas[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const Vector2D<float> *getThetas(int sourceNumber) const 						{ return &(m_thetas[sourceNumber][0]); }
	const Vector2D<float> *getThetas(int sourceNumber, int imageNumber) const				{ return &(m_thetas[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getDerivativesXX(int sourceNumber) const							{ return &(m_axx[sourceNumber][0]); }
	const float *getDerivativesXX(int sourceNumber, int imageNumber) const					{ return &(m_axx[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getDerivativesYY(int sourceNumber) const							{ return &(m_ayy[sourceNumber][0]); }
	const float *getDerivativesYY(int sourceNumber, int imageNumber) const					{ return &(m_ayy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getDerivativesXY(int sourceNumber) const							{ return &(m_axy[sourceNumber][0]); }
	const float *getDerivativesXY(int sourceNumber, int imageNumber) const					{ return &(m_axy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getInverseMagnifications(int sourceNumber) const						{ return &(m_inverseMagnifications[sourceNumber][0]); }
	const float *getInverseMagnifications(int sourceNumber, int imageNumber) const				{ return &(m_inverseMagnifications[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getShearComponents1(int sourceNumber) const						{ return &(m_shearComponent1[sourceNumber][0]); }
	const float *getShearComponents1(int sourceNumber, int imageNumber) const				{ return &(m_shearComponent1[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getShearComponents2(int sourceNumber) const						{ return &(m_axy[sourceNumber][0]); }
	const float *getShearComponents2(int sourceNumber, int imageNumber) const				{ return &(m_axy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getConvergence(int sourceNumber) const							{ return &(m_convergence[sourceNumber][0]); }
	const float *getConvergence(int sourceNumber, int imageNumber) const					{ return &(m_convergence[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	float getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const;

	const float *getLensPotential(int sourceNumber) const { return &(m_potentials[sourceNumber][0]); }
	const float *getLensPotential(int sourceNumber, int imageNumber) const { return &(m_potentials[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
private:
	DeflectionMatrix *m_pDeflectionMatrix;
	bool m_init, m_initializing;

	std::vector<std::vector<int> > m_deflectionIndices;
	std::vector<std::vector<int> > m_derivativeIndices;
	std::vector<std::vector<int> > m_potentialIndices;

	std::vector<std::vector<Vector2D<double> > > m_originalPoints; 

	bool m_useBaseLens;
	std::vector<std::vector<float> > m_baseAxx, m_baseAyy, m_baseAxy;
	std::vector<std::vector<float> > m_basePotentials;
	std::vector<std::vector<Vector2D<float> > > m_baseAlphas;
	std::vector<std::vector<Vector2D<double> > > m_baseAlphasUnscaled;
	std::vector<std::vector<double> > m_basePotentialsUnscaled;

	double m_Dd, m_zd;
	double m_angularScale;
	double m_massSheetScale;
	float m_timeDelayScale;

	std::vector<std::vector<Vector2D<float> > > m_thetas;
	std::vector<std::vector<Vector2D<float> > > m_subDeflectionAngles;
	std::vector<std::vector<float> > m_subDeflectionDerivatives[3];
	std::vector<std::vector<float> > m_subPotentialValues;

	std::vector<std::vector<Vector2D<float> > > m_betas;
	std::vector<std::vector<Vector2D<float> > > m_alphas;
	std::vector<std::vector<float> > m_axx, m_ayy, m_axy;
	std::vector<std::vector<float> > m_potentials;
	std::vector<std::vector<float> > m_inverseMagnifications;
	std::vector<std::vector<float> > m_shearComponent1;
	std::vector<std::vector<float> > m_convergence;

	std::vector<std::vector<Vector2D<float> > > m_sheetAlphas;
	std::vector<std::vector<float> > m_sheetAxx, m_sheetAyy, m_sheetAxy;
	std::vector<std::vector<float> > m_sheetPotentials;
	std::vector<std::vector<Vector2D<double> > > m_sheetAlphasUnscaled;
	std::vector<std::vector<double> > m_sheetPotentialsUnscaled;
	bool m_useMassSheet;

	std::vector<float> m_tmpBuffer;
	std::vector<float> m_oneVector;
	std::vector<bool> m_trueFlags;
};

inline float BackProjectMatrix::getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const
{
	Vector2D<float> diff = beta - m_thetas[sourceNumber][m_offsets[sourceNumber][imageNumber]+pointNumber];
	
	return m_timeDelayScale*(0.5f*diff.getLengthSquared() - m_potentials[sourceNumber][m_offsets[sourceNumber][imageNumber]+pointNumber])/m_distanceFractions[sourceNumber];
}

} // end namespace
