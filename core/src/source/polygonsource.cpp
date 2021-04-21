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
#include "polygonsource.h"

using namespace std;

namespace grale
{
	
PolygonSource::PolygonSource(Vector2D<double> angularpos, const Polygon2D<double> &p, double rotangle, double brightnessScale) : SourceImage(SourceImage::Polygon, angularpos, rotangle, brightnessScale)
{
	polygon = p;
	calcMaxRadius();
}

PolygonSource::~PolygonSource()
{
}

unique_ptr<SourceImage> PolygonSource::createCopy() const
{
	return make_unique<PolygonSource>(getAngularPosition(), polygon, getAngle(), getBrightnessScale());
}
	
double PolygonSource::getIntensityInternal(Vector2Dd diff) const
{
	if (polygon.isInside(diff))
		return 1.0;
	return 0;
}

void PolygonSource::calcMaxRadius()
{
	const Vector2Dd *pPoints = polygon.getPoints();
	int numPoints = polygon.getNumberOfPoints();

	m_maxRadius = pPoints[0].getLength();
	for (int i = 1 ; i < numPoints ; i++)
		m_maxRadius = MAX(m_maxRadius, pPoints[i].getLength());
}

} // end namespace

