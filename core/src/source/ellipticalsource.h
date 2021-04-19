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

#ifndef GRALE_ELLIPTICALSOURCE_H

#define GRALE_ELLIPTICALSOURCE_H

#include "graleconfig.h"
#include "sourceimage.h"
#include "constants.h"

namespace grale
{

class GRALE_IMPORTEXPORT EllipticalSource : public SourceImage
{
public:
	EllipticalSource(Vector2Dd angularpos, double angular_axis, double excentricity, double rotangle, double brightnessScale);
	~EllipticalSource();
	void setFade(bool f)							{ fade = f; }
    bool getFade() const                            { return fade; }
	std::unique_ptr<SourceImage> createCopy() const override;
protected:
	double getIntensityInternal(Vector2Dd diff) const;
	double getMaxRadius() const						{ return A; }
private:
	double A,B,e;
	bool fade;
};

inline double EllipticalSource::getIntensityInternal(Vector2Dd diff) const
{
	// rescale 
	Vector2Dd diff3(diff.getX()/A,diff.getY()/B);
	
	if (diff3.getLengthSquared() <= 1.0)
	{
		double val = 1.0;
		
		if (fade)
			val *= (1.0L-(diff3.getLengthSquared()));
		
		return val;
	}
	return 0;
}

} // end namespace

#endif // GRALE_ELLIPTICALSOURCE_H

