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

#ifndef GRALE_IMAGESBACKPROJECTOR_H

#define GRALE_IMAGESBACKPROJECTOR_H

/**
 * \file imagesbackprojector.h
 */

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include "vector2d.h"
#include <errut/errorbase.h>
#include <vector>
#include <list>
#include <string>
#include <memory>

namespace grale
{

class ImagesDataExtended;
class GravitationalLens;
class LensFitnessObject;

/** Implements the ProjectedImagesInterface interface, and takes a GravitationalLens based
 *  lens model as input. */
class GRALE_IMPORTEXPORT ImagesBackProjector : public ProjectedImagesInterface, public errut::ErrorBase
{
public:
	/** Contructor of the class.
	 *  Contructor of the class.
	 *  \param lens The lens model
	 *  \param images The list of images data (actual lensed images, null space etc)
	 *  \param z_d The redshift of the lens, important if time delays will be calculated.
	 *  \param copyLens If true, a copy of the lens model is created and stored internally.
	 *                  Otherwise, a pointer to the lens model is stored.
	 */
	ImagesBackProjector(const std::shared_ptr<GravitationalLens> &lens, const std::list<ImagesDataExtended *> &images, double z_d); 
	~ImagesBackProjector();
	
	double getLensDistance() const;
	double getLensRedshift() const							{ return m_zd; }
	double getAngularScale() const							{ return m_angularScale; }
	const Vector2D<float> *getBetas(int sourceNumber) const				{ checkBetas(sourceNumber); return &(m_betas[sourceNumber][0]); }
	const Vector2D<float> *getBetas(int sourceNumber, int imageNumber) const	{ checkBetas(sourceNumber); return &(m_betas[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const Vector2D<float> *getAlphas(int sourceNumber) const				{ checkAlphas(sourceNumber); return &(m_alphas[sourceNumber][0]); }
	const Vector2D<float> *getAlphas(int sourceNumber, int imageNumber) const	{ checkAlphas(sourceNumber); return &(m_alphas[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const Vector2D<float> *getThetas(int sourceNumber) const			{ return &(m_thetas[sourceNumber][0]); }
	const Vector2D<float> *getThetas(int sourceNumber, int imageNumber) const	{ return &(m_thetas[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getDerivativesXX(int sourceNumber) const				{ checkDerivatives(sourceNumber); return &(m_axx[sourceNumber][0]); }
	const float *getDerivativesXX(int sourceNumber, int imageNumber) const		{ checkDerivatives(sourceNumber); return &(m_axx[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getDerivativesYY(int sourceNumber) const 				{ checkDerivatives(sourceNumber); return &(m_ayy[sourceNumber][0]); }
	const float *getDerivativesYY(int sourceNumber, int imageNumber) const 		{ checkDerivatives(sourceNumber); return &(m_ayy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getDerivativesXY(int sourceNumber) const 				{ checkDerivatives(sourceNumber); return &(m_axy[sourceNumber][0]); }
	const float *getDerivativesXY(int sourceNumber, int imageNumber) const		{ checkDerivatives(sourceNumber); return &(m_axy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }

	const float *getSecondDerivativesXXX(int sourceNumber) const						{ checkSecondDerivatives(sourceNumber); return &(m_axxx[sourceNumber][0]); }
	const float *getSecondDerivativesXXX(int sourceNumber, int imageNumber) const		{ checkSecondDerivatives(sourceNumber); return &(m_axxx[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getSecondDerivativesYYY(int sourceNumber) const						{ checkSecondDerivatives(sourceNumber); return &(m_ayyy[sourceNumber][0]); }
	const float *getSecondDerivativesYYY(int sourceNumber, int imageNumber) const		{ checkSecondDerivatives(sourceNumber); return &(m_ayyy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getSecondDerivativesXXY(int sourceNumber) const						{ checkSecondDerivatives(sourceNumber); return &(m_axxy[sourceNumber][0]); }
	const float *getSecondDerivativesXXY(int sourceNumber, int imageNumber) const		{ checkSecondDerivatives(sourceNumber); return &(m_axxy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getSecondDerivativesYYX(int sourceNumber) const						{ checkSecondDerivatives(sourceNumber); return &(m_ayyx[sourceNumber][0]); }
	const float *getSecondDerivativesYYX(int sourceNumber, int imageNumber) const		{ checkSecondDerivatives(sourceNumber); return &(m_ayyx[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }

	const float *getInverseMagnifications(int sourceNumber) const			{ checkInvMag(sourceNumber); return &(m_invMag[sourceNumber][0]); }
	const float *getInverseMagnifications(int sourceNumber, int imageNumber) const	{ checkInvMag(sourceNumber); return &(m_invMag[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getShearComponents1(int sourceNumber) const 			{ checkShear(sourceNumber); return &(m_shearComponents1[sourceNumber][0]); }
	const float *getShearComponents1(int sourceNumber, int imageNumber) const	{ checkShear(sourceNumber); return &(m_shearComponents1[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getShearComponents2(int sourceNumber) const			{ checkDerivatives(sourceNumber); return &(m_axy[sourceNumber][0]); }
	const float *getShearComponents2(int sourceNumber, int imageNumber) const 	{ checkDerivatives(sourceNumber); return &(m_axy[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	const float *getConvergence(int sourceNumber) const				{ checkConvergence(sourceNumber); return &(m_convergence[sourceNumber][0]); }
	const float *getConvergence(int sourceNumber, int imageNumber) const		{ checkConvergence(sourceNumber); return &(m_convergence[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }

	const float *getLensPotential(int sourceNumber) const { checkPotential(sourceNumber); return &(m_potential[sourceNumber][0]); }
	const float *getLensPotential(int sourceNumber, int imageNumber) const { checkPotential(sourceNumber); return &(m_potential[sourceNumber][m_offsets[sourceNumber][imageNumber]]); }
	float getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const;
private:
	void checkBetas(int sourceNumber) const;
	void checkAlphas(int sourceNumber) const;
	void checkDerivatives(int sourceNumber) const;
	void checkSecondDerivatives(int sourceNumber) const;
	void checkInvMag(int sourceNumber) const;
	void checkShear(int sourceNumber) const;
	void checkConvergence(int sourceNumber) const;
	void checkPotential(int sourceNumber) const;

	mutable std::vector<std::vector<Vector2D<float> > > m_betas;
	mutable std::vector<std::vector<Vector2D<float> > > m_alphas;
	std::vector<std::vector<Vector2D<float> > > m_thetas;
	std::vector<std::vector<Vector2D<double> > > m_originalThetas;
	mutable std::vector<std::vector<float> > m_axx, m_ayy, m_axy;
	mutable std::vector<std::vector<float> > m_axxx, m_ayyy, m_axxy, m_ayyx;
	mutable std::vector<std::vector<float> > m_invMag;
	mutable std::vector<std::vector<float> > m_shearComponents1;
	mutable std::vector<std::vector<float> > m_convergence;
	mutable std::vector<std::vector<float> > m_potential;
	double m_angularScale;

	std::shared_ptr<GravitationalLens> m_pLens;
	double m_zd;
};

} // end namespace

#endif // GRALE_IMAGESBACKPROJECTOR_H

