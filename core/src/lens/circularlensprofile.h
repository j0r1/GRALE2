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

#ifndef GRALE_CIRCULARLENSPROFILE_H

#define GRALE_CIRCULARLENSPROFILE_H

namespace grale
{

class GRALE_IMPORTEXPORT CircularLensProfile
{
public:
	CircularLensProfile()										{ }
	virtual ~CircularLensProfile()									{ }
	virtual double getMassInside(double theta) const = 0;
	virtual double getSurfaceMassDensity(double theta) const = 0;
	virtual double getSurfaceMassDensityDerivativeOverTheta(double theta) const = 0;
};

} // end namespace

#endif // GRALE_CIRCULARLENSPROFILE_H
