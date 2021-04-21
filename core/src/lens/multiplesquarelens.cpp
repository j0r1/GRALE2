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
#include "multiplesquarelens.h"
#include "constants.h"
#include <iostream>

#include <stdio.h>

namespace grale
{

std::unique_ptr<GravitationalLensParams> MultipleSquareLensParams::createCopy() const
{
	return std::make_unique<MultipleSquareLensParams>(m_lensInfo);
}

bool MultipleSquareLensParams::write(serut::SerializationInterface &si) const
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
			setErrorString(std::string("Error writing square lens parameters: ") + si.getErrorString());
			return false;
		}
	}
	return true;
}

bool MultipleSquareLensParams::read(serut::SerializationInterface &si)
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
			setErrorString(std::string("Error reading square lens parameters: ") + si.getErrorString());
			return false;
		}
		
		Vector2D<double> v(array[2],array[3]);
		SquareLensInfo p(array[0],array[1],v);
		m_lensInfo.push_back(p);
	}
	
	return true;
}

MultipleSquareLens::MultipleSquareLens() : GravitationalLens(GravitationalLens::MultipleSquares)
{
}

MultipleSquareLens::~MultipleSquareLens()
{
}

bool MultipleSquareLens::processParameters(const GravitationalLensParams *params)
{
	const MultipleSquareLensParams *p = dynamic_cast<const MultipleSquareLensParams *>(params);
	if (!p)
	{
		setErrorString("Parameters are not of type 'MultipleSquareLensParams'");
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
		lensinfo[i] = SquareLensInfo((*it).getMass()/totalmass,(*it).getAngularWidth()/scalefactor,(*it).getAngularPosition());
	
	return true;
}

bool MultipleSquareLens::getAlphaVector(Vector2D<double> theta,Vector2D<double> *alpha) const
{	
	Vector2D<double> a(0.0,0.0);
	
	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> scaledtheta = diff/scalefactor;
		double m = lensinfo[i].getMass();
		double w = lensinfo[i].getAngularWidth();
		double w2 = w*w;
		double rawd2 = w/2.0;

		double ax = deflectionFunction1(scaledtheta.getX()-rawd2,scaledtheta.getY()-rawd2)
			  - deflectionFunction1(scaledtheta.getX()+rawd2,scaledtheta.getY()-rawd2)
			  - deflectionFunction1(scaledtheta.getX()-rawd2,scaledtheta.getY()+rawd2)
			  + deflectionFunction1(scaledtheta.getX()+rawd2,scaledtheta.getY()+rawd2);
		double ay = deflectionFunction2(scaledtheta.getX()-rawd2,scaledtheta.getY()-rawd2)
			  - deflectionFunction2(scaledtheta.getX()+rawd2,scaledtheta.getY()-rawd2)
			  - deflectionFunction2(scaledtheta.getX()-rawd2,scaledtheta.getY()+rawd2)
			  + deflectionFunction2(scaledtheta.getX()+rawd2,scaledtheta.getY()+rawd2);

		a += Vector2D<double>(ax,ay)*m/w2;
	}
	
	a *= scalefactor;
	*alpha = a;
	return true;
}

double MultipleSquareLens::getSurfaceMassDensity(Vector2D<double> theta) const
{
	double sum = 0;

	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> scaledtheta = diff/scalefactor;
		double w = lensinfo[i].getAngularWidth();

		if (scaledtheta.getX() > -w/2.0 && scaledtheta.getX() < w/2.0 &&
		    scaledtheta.getY() > -w/2.0 && scaledtheta.getY() < w/2.0)
		{
			
			double m = lensinfo[i].getMass();
			double w2 = w*w;
			double dens = m/w2;

			sum += dens;
		}
	}
	return sum*totalmass/(getLensDistance()*getLensDistance()*scalefactor*scalefactor);
}

bool MultipleSquareLens::getAlphaVectorDerivatives(Vector2D<double> theta, double &axx, double &ayy, double &axy) const
{
	double a11 = 0;
	double a22 = 0;
	double a12 = 0;

	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> t = diff/scalefactor;
		double m = lensinfo[i].getMass();
		double w = lensinfo[i].getAngularWidth();
		double w2 = w*w;
		double rawd2 = w/2.0;
		
		a11 += m*(1.0/w2)*(deriv11(t.getX()-rawd2,t.getY()-rawd2)
							 -deriv11(t.getX()+rawd2,t.getY()-rawd2)
							 -deriv11(t.getX()-rawd2,t.getY()+rawd2)
							 +deriv11(t.getX()+rawd2,t.getY()+rawd2));
		a22 += m*(1.0/w2)*(deriv22(t.getX()-rawd2,t.getY()-rawd2)
							 -deriv22(t.getX()+rawd2,t.getY()-rawd2)
							 -deriv22(t.getX()-rawd2,t.getY()+rawd2)
							 +deriv22(t.getX()+rawd2,t.getY()+rawd2));
		a12 += m*(1.0/w2)*(deriv12(t.getX()-rawd2,t.getY()-rawd2)
							 -deriv12(t.getX()+rawd2,t.getY()-rawd2)
							 -deriv12(t.getX()-rawd2,t.getY()+rawd2)
							 +deriv12(t.getX()+rawd2,t.getY()+rawd2));
	}
	axx = a11;
	ayy = a22;
	axy = a12;
	return true;
}

bool MultipleSquareLens::getProjectedPotential(double D_s, double D_ds, Vector2D<double> theta, double *potentialval) const
{
	double sum = 0;
	
	for (int i = 0 ; i < numlenses ; i++)
	{
		Vector2D<double> diff = theta-lensinfo[i].getAngularPosition();
		Vector2D<double> t = diff/scalefactor;
		double m = lensinfo[i].getMass();
		double w = lensinfo[i].getAngularWidth();
		double w2 = w*w;
		double rawd2 = w/2.0;
		
		double p = potential(t.getX()-rawd2,t.getY()-rawd2)-potential(t.getX()+rawd2,t.getY()-rawd2)
			 - potential(t.getX()-rawd2,t.getY()+rawd2)+potential(t.getX()+rawd2,t.getY()+rawd2);

		sum += p*m/w2;
	}
	
	*potentialval = (D_ds/D_s)*scalefactor*scalefactor*sum;	
	return true;
}

} // end namespace

