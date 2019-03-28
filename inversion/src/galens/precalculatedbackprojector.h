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

#ifndef GRALE_PRECALCULATEDBACKPROJECTOR_H

#define GRALE_PRECALCULATEDBACKPROJECTOR_H

/**
 * \file precalculatedbackprojector.h
 */

#include "graleconfig.h"
#include "projectedimagesinterface.h"
#include "vector2d.h"
#include <errut/errorbase.h>
#include <assert.h>
#include <vector>
#include <limits>
#include <string>

namespace grale
{

class ImagesData;

class GRALE_IMPORTEXPORT PreCalculatedBackProjector : public ProjectedImagesInterface, public errut::ErrorBase
{
public:
	PreCalculatedBackProjector();
	~PreCalculatedBackProjector();

	bool init(const std::vector<ImagesData *> &images, const std::vector<ImagesData *> &correspondingSources); 
	
	double getLensDistance() const														{ return std::numeric_limits<double>::quiet_NaN(); }
	double getLensRedshift() const														{ return std::numeric_limits<double>::quiet_NaN(); }
	double getAngularScale() const														{ return m_angularScale; }
	const Vector2D<float> *getBetas(int sourceNumber) const;
	const Vector2D<float> *getBetas(int sourceNumber, int imageNumber) const;
	const Vector2D<float> *getThetas(int sourceNumber) const;
	const Vector2D<float> *getThetas(int sourceNumber, int imageNumber) const;
	const Vector2D<float> *getAlphas(int sourceNumber) const							{ return nullptr; }
	const Vector2D<float> *getAlphas(int sourceNumber, int imageNumber) const			{ return nullptr; }
	const float *getDerivativesXX(int sourceNumber) const								{ return nullptr; }
	const float *getDerivativesXX(int sourceNumber, int imageNumber) const				{ return nullptr; }
	const float *getDerivativesYY(int sourceNumber) const								{ return nullptr; }
	const float *getDerivativesYY(int sourceNumber, int imageNumber) const				{ return nullptr; }
	const float *getDerivativesXY(int sourceNumber) const								{ return nullptr; }
	const float *getDerivativesXY(int sourceNumber, int imageNumber) const				{ return nullptr; }
	const float *getInverseMagnifications(int sourceNumber) const						{ return nullptr; }
	const float *getInverseMagnifications(int sourceNumber, int imageNumber) const		{ return nullptr; }
	const float *getShearComponents1(int sourceNumber) const							{ return nullptr; }
	const float *getShearComponents1(int sourceNumber, int imageNumber) const			{ return nullptr; }
	const float *getShearComponents2(int sourceNumber) const							{ return nullptr; }
	const float *getShearComponents2(int sourceNumber, int imageNumber) const			{ return nullptr; }
	const float *getConvergence(int sourceNumber) const									{ return nullptr; }
	const float *getConvergence(int sourceNumber, int imageNumber) const				{ return nullptr; }

	const float *getLensPotential(int sourceNumber) const								{ return nullptr; }
	const float *getLensPotential(int sourceNumber, int imageNumber) const				{ return nullptr; }
	float getTimeDelay(int sourceNumber, int imageNumber, int pointNumber, Vector2D<float> beta) const { return std::numeric_limits<float>::quiet_NaN(); }
private:
	std::vector<std::vector<Vector2D<float> > > m_betas, m_thetas;
	double m_angularScale;
};

inline const Vector2D<float> *PreCalculatedBackProjector::getBetas(int sourceNumber) const
{
	assert(sourceNumber < m_betas.size());
	return &(m_betas[sourceNumber][0]);
}

inline const Vector2D<float> *PreCalculatedBackProjector::getBetas(int sourceNumber, int imageNumber) const
{
	assert(sourceNumber < m_betas.size() && sourceNumber < m_offsets.size());
	assert(imageNumber < m_offsets[sourceNumber].size());

	int offset = m_offsets[sourceNumber][imageNumber];
	assert(offset < m_betas[sourceNumber].size());
	return &(m_betas[sourceNumber][offset]);
}

inline const Vector2D<float> *PreCalculatedBackProjector::getThetas(int sourceNumber) const
{
	assert(sourceNumber < m_thetas.size());
	return &(m_thetas[sourceNumber][0]);
}

inline const Vector2D<float> *PreCalculatedBackProjector::getThetas(int sourceNumber, int imageNumber) const
{
	assert(sourceNumber < m_thetas.size() && sourceNumber < m_offsets.size());
	assert(imageNumber < m_offsets[sourceNumber].size());

	int offset = m_offsets[sourceNumber][imageNumber];
	assert(offset < m_thetas[sourceNumber].size());
	return &(m_thetas[sourceNumber][offset]);
}

} // end namespace

#endif // GRALE_PRECALCULATEDBACKPROJECTOR_H

