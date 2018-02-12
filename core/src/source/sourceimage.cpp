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
#include "sourceimage.h"

#include "debugnew.h"

namespace grale
{

SourceImage::SourceImage(SourceType t,Vector2D<double> pos, double angle, double brightnessScale)
{
	SourceImage::pos = pos;
	SourceImage::angle = angle;
	m_brightnessScale = brightnessScale;
	stype = t;
	theta = (angle/180.0)*CONST_PI;
}

SourceImage::~SourceImage()
{
}

} // end namespace

