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

#ifndef GRALE_CIRCULARSOURCE_H

#define GRALE_CIRCULARSOURCE_H

#include "graleconfig.h"
#include "sourceimage.h"

namespace grale
{

class GRALE_IMPORTEXPORT CircularSource : public SourceImage
{
public:
	CircularSource(Vector2Dd angularpos, double angularradius, double brightnessScale);
	~CircularSource();
	SourceImage *createCopy() const;
	double getAngularRadius() const						{ return radius; }
	void setAngularRadius(double r)						{ radius = r; radius2 = r*r; }
	void setFade(bool f)								{ fade = f; }
	bool getFade() const								{ return fade; }
protected:
	double getIntensityInternal(Vector2Dd diff) const;
	double getMaxRadius() const						{ return radius; }
private:
	double radius,radius2;
	bool fade;
};

inline double CircularSource::getIntensityInternal(Vector2Dd diff) const
{
	if (diff.getLengthSquared() <= radius2)
	{
		double val = 1.0;
		
		if (fade)
			val *= (1.0-(diff.getLengthSquared()/radius2));
		
		return val;
	}
	return 0;
}

} // end namespace

#endif // GRALE_CIRCULARSOURCE_H

