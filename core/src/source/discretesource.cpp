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
#include "discretesource.h"
#include <string.h>

#include "debugnew.h"

using namespace std;

namespace grale
{

DiscreteSource::DiscreteSource(Vector2D<double> angularpos, double angle, double brightnessScale) : SourceImage(SourceImage::Discrete, angularpos, angle, brightnessScale)
{
	m_maxRadius = 0;
	m_numX = 0;
	m_numY = 0;
	m_width = 0;
	m_height = 0;
}

DiscreteSource::~DiscreteSource()
{
}

bool DiscreteSource::setData(const std::vector<double> &data, int numX, int numY, double angularWidth, double angularHeight)
{
	if (numX < 1 || numY < 1)
	{
		setErrorString("Number of pixels in X and Y directions must be at least one");
		return false;
	}
	if (numX*numY != (int)data.size())
	{
		setErrorString("Total number of pixels doesn't match data length");
		return false;
	}
	if (angularWidth <= 0 || angularHeight <= 0)
	{
		setErrorString("Angular width and height of image must be positive");
		return false;
	}

	m_data = data;
	m_numX = numX;
	m_numY = numY;
	m_width = angularWidth;
	m_height = angularHeight;

	calcMaxRadius();
	return true;
}

unique_ptr<SourceImage> DiscreteSource::createCopy() const
{
	auto pNewSrc = make_unique<DiscreteSource>(getAngularPosition(), getAngle(), getBrightnessScale());

	pNewSrc->m_data = m_data;
	pNewSrc->m_numX = m_numX;
	pNewSrc->m_numY = m_numY;
	pNewSrc->m_width = m_width;
	pNewSrc->m_height = m_height;
	pNewSrc->calcMaxRadius();
	return pNewSrc;
}

double DiscreteSource::getIntensityInternal(Vector2D<double> diff) const
{
	double x = diff.getX()+m_width/2.0;
	double y = diff.getY()+m_height/2.0;

	if (x < 0 || y < 0 || x > m_width || y > m_height)
		return 0;

	int xPix = (int)((x/m_width)*m_numX);
	int yPix = m_numY-1-(int)((y/m_height)*m_numY);
	if (xPix < 0 || xPix >= m_numX || yPix < 0 || yPix >= m_numY)
		return 0;
	
	int pos = xPix + yPix*m_numX;
	assert(pos >= 0 && pos < (int)m_data.size());
	return m_data[pos];
}

void DiscreteSource::calcMaxRadius()
{
	double x1 = -m_width/2.0;
	double y1 = -m_height/2.0;
	double x2 = m_width/2.0;
	double y2 = m_height/2.0;

	Vector2D<double> c1(x1,y1);
	Vector2D<double> c2(x1,y2);
	Vector2D<double> c3(x2,y1);
	Vector2D<double> c4(x2,y2);

	m_maxRadius = c1.getLength();
	m_maxRadius = MAX(m_maxRadius, c2.getLength());
	m_maxRadius = MAX(m_maxRadius, c3.getLength());
	m_maxRadius = MAX(m_maxRadius, c4.getLength());
}

} // end namespace


