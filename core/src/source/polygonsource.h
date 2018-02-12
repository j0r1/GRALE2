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

#ifndef GRALE_POLYGONSOURCE_H

#define GRALE_POLYGONSOURCE_H

#include "graleconfig.h"
#include "sourceimage.h"
#include "polygon2d.h"

namespace grale
{

class GRALE_IMPORTEXPORT PolygonSource : public SourceImage
{
public:
	PolygonSource(Vector2Dd angularpos, const Polygon2D<double> &p, double rotangle, double brightnessScale);
	~PolygonSource();

	SourceImage *createCopy() const;
protected:
	double getIntensityInternal(Vector2Dd diff) const;
	double getMaxRadius() const						{ return m_maxRadius; }
private:
	void calcMaxRadius();

	Polygon2D<double> polygon;
	double m_maxRadius;
};
	
} // end namespace

#endif // GRALE_POLYGONSOURCE_H

