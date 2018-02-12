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

#ifndef GRALE_CONSTANTS_H

#define GRALE_CONSTANTS_H

namespace grale
{

const double MASS_SOLAR = 1.98855e30;
const double DIST_PC = 3.0856775714409184e16;
const double DIST_KPC = 1000.0*DIST_PC;
const double DIST_MPC = 1000.0*DIST_KPC;
const double DIST_GPC = 1000.0*DIST_MPC;
const double SPEED_C = 299792458.0;
const double CONST_G = 6.67300e-11;
const double CONST_PI = 3.14159265358979323846264338327950288419;
const double ANGLE_DEGREE = CONST_PI/180.0;
const double ANGLE_ARCMIN = ANGLE_DEGREE/60.0;
const double ANGLE_ARCSEC = ANGLE_ARCMIN/60.0;
const double ANGLE_MILLIARCSEC = ANGLE_ARCSEC/1000.0;
const double ANGLE_MICROARCSEC = ANGLE_MILLIARCSEC/1000.0;

}

#endif // GRALE_CONSTANTS_H

