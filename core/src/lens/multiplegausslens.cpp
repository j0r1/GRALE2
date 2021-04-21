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
#include "multiplegausslens.h"
#include "constants.h"
#include <iostream>

#include "debugnew.h"

namespace grale
{

std::unique_ptr<GravitationalLensParams> MultipleGaussLensParams::createCopy() const
{
	return std::make_unique<MultipleGaussLensParams>(m_lensInfo);
}

bool MultipleGaussLensParams::write(serut::SerializationInterface &si) const
{
	int32_t numlenses = (int32_t)m_lensInfo.size();
	
	if (!si.writeInt32(numlenses))
	{
		setErrorString(std::string("Error writing number of lenses: ") + si.getErrorString());
		return false;
	}

	for (auto it = m_lensInfo.begin() ; it != m_lensInfo.end() ; it++)
	{
		double array[4];
		
		array[0] = (*it).getMass();
		array[1] = (*it).getAngularWidth();
		array[2] = (*it).getAngularPosition().getX();
		array[3] = (*it).getAngularPosition().getY();
		
		if (!si.writeDoubles(array,4))
		{
			setErrorString(std::string("Error writing gauss lens parameters: ") + si.getErrorString());
			return false;
		}
	}
	return true;
}

bool MultipleGaussLensParams::read(serut::SerializationInterface &si)
{
	int32_t numlenses,i;

	if (!si.readInt32(&numlenses))
	{
		setErrorString(std::string("Error reading number of lenses: ") + si.getErrorString());
		return false;
	}
	
	m_lensInfo.clear();
	for (i = 0 ; i < numlenses ; i++)
	{
		double array[4];

		if (!si.readDoubles(array,4))
		{
			setErrorString(std::string("Error reading gauss lens parameters: ") + si.getErrorString());
			return false;
		}
		
		Vector2D<double> v(array[2],array[3]);
		GaussLensInfo p(array[0],array[1],v);
		m_lensInfo.push_back(p);
	}
	
	return true;
}

MultipleGaussLens::MultipleGaussLens() : GravitationalLens(GravitationalLens::MultipleGaussians)
{
}

MultipleGaussLens::~MultipleGaussLens()
{
}

bool MultipleGaussLens::processParameters(const GravitationalLensParams *params)
{
	const MultipleGaussLensParams *p = dynamic_cast<const MultipleGaussLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'MultipleGaussLensParams'");
		return false;
	}
	
	int i;

	numlenses = p->getLensInfo().size();
	lensinfo.resize(numlenses);

	totalmass = 0;
	auto it = p->getLensInfo().begin();
	for (i = 0 ; i < numlenses ; i++,it++)
		totalmass += (*it).getMass();

	totalmass = ABS(totalmass);
	scalefactor = SQRT((4.0*CONST_G*totalmass)/(SPEED_C*SPEED_C*getLensDistance()));
		
	for (i = 0, it = p->getLensInfo().begin() ; i < numlenses ; i++,it++)
		lensinfo[i] = GaussLensInfo((*it).getMass()/totalmass,(*it).getAngularWidth()/scalefactor,(*it).getAngularPosition());
	
	return true;
}

bool MultipleGaussLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *alpha) const
{	
	Vector2D<double> sum(0,0);
	
	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[i].getAngularWidth();
		double m = lensinfo[i].getMass();
		double l2 = scaleddiff.getLengthSquared();

		sum += m*((1.0-EXP(-l2/(2.0*w*w)))/l2)*scaleddiff;
	}

	*alpha = sum*scalefactor;
	return true;
}

double MultipleGaussLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double sum = 0;
	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[i].getAngularWidth();
		double m = lensinfo[i].getMass();

		sum += m*(1.0/(w*w))*EXP(-scaleddiff.getLengthSquared()/(2.0*w*w));;
	}
	return (sum*SPEED_C*SPEED_C)/(8.0*CONST_PI*CONST_G*getLensDistance());
}

double MultipleGaussLens::getDeriv11(Vector2D<double> theta) const
{
	double deriv = 0;
	
	for (int k = 0 ; k < numlenses ; k++)
	{
		Vector2D<double> diff = theta-lensinfo[k].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[k].getAngularWidth();
		double m = lensinfo[k].getMass();
		double x = scaleddiff.getX();
		double y = scaleddiff.getY();
		double e = EXP(-scaleddiff.getLengthSquared()/(2.0*w*w));
		
		deriv += m/(scaleddiff.getLengthSquared())*((x/w)*(x/w)*e+(1.0-e)*(y*y-x*x)/scaleddiff.getLengthSquared());
	}
	return deriv;
}

double MultipleGaussLens::getDeriv22(Vector2D<double> theta) const
{
	double deriv = 0;
	
	for (int k = 0 ; k < numlenses ; k++)
	{
		Vector2D<double> diff = theta-lensinfo[k].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[k].getAngularWidth();
		double m = lensinfo[k].getMass();
		double x = scaleddiff.getX();
		double y = scaleddiff.getY();
		double e = EXP(-scaleddiff.getLengthSquared()/(2.0*w*w));
		
		deriv += m/(scaleddiff.getLengthSquared())*((y/w)*(y/w)*e+(1.0-e)*(x*x-y*y)/scaleddiff.getLengthSquared());
	}
	return deriv;
}

double MultipleGaussLens::getDeriv12(Vector2D<double> theta) const
{
	double deriv = 0;
	
	for (int k = 0 ; k < numlenses ; k++)
	{
		Vector2D<double> diff = theta-lensinfo[k].getAngularPosition();
		Vector2D<double> scaleddiff = diff/scalefactor;
		double w = lensinfo[k].getAngularWidth();
		double m = lensinfo[k].getMass();
		double x = scaleddiff.getX();
		double y = scaleddiff.getY();
		double e = EXP(-scaleddiff.getLengthSquared()/(2.0*w*w));
		
		deriv += m/(scaleddiff.getLengthSquared())*((x/w)*(y/w)*e-2*(1.0-e)*(y*x)/scaleddiff.getLengthSquared());
	}

	return deriv;
}

bool MultipleGaussLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *potentialval) const
{
	// TODO: implement this
	return false;
}

bool MultipleGaussLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	axx = getDeriv11(theta);
	ayy = getDeriv22(theta);
	axy = getDeriv12(theta);
	return true;
}

} // end namespace

