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
#include "profilelens.h"
#include "constants.h"

#include "debugnew.h"

namespace grale
{
	
GravitationalLensParams *ProfileLensParams::createCopy() const
{ 
	return new ProfileLensParams(m_endRadius, m_profile); 
}

bool ProfileLensParams::write(serut::SerializationInterface &si) const
{ 
	if (!si.writeInt32(m_profile.size()))
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	if (!si.writeDoubles(m_profile)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	if (!si.writeDouble(m_endRadius)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

bool ProfileLensParams::read(serut::SerializationInterface &si)
{ 
	int32_t numPoints; 
	
	if (!si.readInt32(&numPoints)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	m_profile.resize(numPoints); 
	if (!si.readDoubles(m_profile)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	if (!si.readDouble(&m_endRadius)) 
	{ 
		setErrorString(si.getErrorString()); 
		return false; 
	} 
	
	return true; 
}

ProfileLens::ProfileLens() : SymmetricLens(GravitationalLens::Profile)
{
}

ProfileLens::~ProfileLens()
{
}

bool ProfileLens::processParameters(const GravitationalLensParams *pLensParams)
{
	const ProfileLensParams *pParams = dynamic_cast<const ProfileLensParams *>(pLensParams);

	if (!pParams)
	{
		setErrorString("Parameters are not of type 'ProfileLensParams'");
		return false;
	}

	m_profile = pParams->getProfile();

	if (m_profile.size() < 2)
	{
		setErrorString("The profile should contain at least two points");
		return false;
	}
	
	m_endRadius = pParams->getEndRadius();
	m_stepSize = m_endRadius/(m_profile.size()-1);

	double factor = 2.0*CONST_PI*getLensDistance()*getLensDistance();
	double sum = 0;
	double theta = m_stepSize;
	double prevValue = 0;
	double potential = 0;
	double prevMassOverTheta = 0;
	
	m_mass.resize(m_profile.size());	
	m_potential.resize(m_profile.size());
	
	m_mass[0] = 0;
	m_potential[0] = 0;
	
	for (int i = 1 ; i < m_mass.size() ; i++)
	{
		double dens = m_profile[i];
		double curValue = dens*theta;

		sum += m_stepSize*(curValue+prevValue)/2.0;

		m_mass[i] = factor*sum;

		potential += m_stepSize*(m_mass[i]/theta+prevMassOverTheta)/2.0;

		m_potential[i] = potential;

		prevValue = curValue;
		prevMassOverTheta = m_mass[i]/theta;
		theta += m_stepSize;
	}

	return true;
}

double ProfileLens::getMassInside(double thetaLength) const
{
	double ratio = thetaLength/m_stepSize;
	int index = (int)ratio;
	
	if (index >= m_mass.size()-2)
		return m_mass[m_mass.size()-1];

	double prevVal = m_mass[index];
	double nextVal = m_mass[index+1];
	double frac = ratio-(double)index;
	
	return (1.0-frac)*prevVal+frac*nextVal;
}

double ProfileLens::getProfileSurfaceMassDensity(double thetaLength) const
{
	double ratio = thetaLength/m_stepSize;
	int index = (int)ratio;
	
	if (index >= m_profile.size()-2)
		return m_profile[m_profile.size()-1];

	double prevVal = m_profile[index];
	double nextVal = m_profile[index+1];
	double frac = ratio-(double)index;
	
	return (1.0-frac)*prevVal+frac*nextVal;
}

bool ProfileLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, 
                                        double *pPotentialValue) const
{
	double thetaLength = theta.getLength();
	double ratio = thetaLength/m_stepSize;
	int index = (int)ratio;
	
	if (index >= m_potential.size()-2)
		return m_potential[m_potential.size()-1];

	double prevVal = m_potential[index];
	double nextVal = m_potential[index+1];
	double frac = ratio-(double)index;
	double factor = ((D_ds/D_s)*(4.0*CONST_G/(SPEED_C*SPEED_C))/getLensDistance());
	
	*pPotentialValue = factor*((1.0-frac)*prevVal+frac*nextVal);
	return true;
}

} // end namespace

