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

#ifndef GRALE_SOURCEIMAGE_H

#define GRALE_SOURCEIMAGE_H

#include "graleconfig.h"
#include "constants.h"
#include "vector2d.h"
#include <errut/errorbase.h>
#include <string>

namespace grale
{

class GRALE_IMPORTEXPORT SourceImage : public errut::ErrorBase
{
public:
	enum SourceType { Circle, Ellipse, Polygon, Discrete, Point };
protected:	
	// angle is specified in degrees
	SourceImage(SourceType t, Vector2Dd angularpos, double angle, double brightnessScale);
public:
	virtual ~SourceImage();

	SourceType getSourceType() const					{ return stype; }
	
	// creates an unattached copy
	virtual SourceImage *createCopy() const = 0;
		
	double getIntensity(Vector2Dd beta) const;
	bool isSourceInRange(Vector2Dd beta, double radius) const; 

	Vector2Dd getAngularPosition() const				{ return pos; }
	void addToAngularPosition(Vector2Dd p)				{ pos += p; }
	void setAngularPosition(Vector2Dd p)				{ pos = p; }
	void addToAngle(double ang)						{ angle += ang; theta = (angle/180.0)*CONST_PI; }
	void setAngle(double ang)						{ angle = ang; theta = (angle/180.0)*CONST_PI; }
	double getAngle() const							{ return angle; }

	virtual double getMaxRadius() const						{ return 0; }
protected:
	double getBrightnessScale() const					{ return m_brightnessScale; }
	
	virtual double getIntensityInternal(Vector2Dd diff) const		{ return 0; }
private:
	Vector2Dd pos;
	double angle,theta;
	SourceType stype;
	double m_brightnessScale;
	
	friend class SourcePlane;
};

inline double SourceImage::getIntensity(Vector2Dd beta) const
{
	Vector2Dd diff = beta - pos;
	Vector2Dd diff2(diff.getX()*COS(theta) + diff.getY()*SIN(theta), 
			      -diff.getX()*SIN(theta) + diff.getY()*COS(theta));
	
	return m_brightnessScale*getIntensityInternal(diff2);
}

inline bool SourceImage::isSourceInRange(Vector2Dd beta, double radius) const
{
	Vector2Dd diff = beta - pos;
	Vector2Dd diff2(diff.getX()*COS(theta) + diff.getY()*SIN(theta), 
			      -diff.getX()*SIN(theta) + diff.getY()*COS(theta));
	
	double d = diff2.getLength() - radius;

	if (d < getMaxRadius())
		return true;
	return false;
}

} // end namespace

#endif // GRALE_SOURCEIMAGE_H

