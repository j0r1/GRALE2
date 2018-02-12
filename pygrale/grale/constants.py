"""This module defines several constants that can be used to express
necessary quantities in.

Throughout the program, basic SI units are assumed, but with the constants
defined in this module it is easy to convert certain quantities to these
units. These constants are:

 - MASS_SUN
 - DIST_PC
 - DIST_KPC
 - DIST_MPC
 - DIST_GPC
 - SPEED_C
 - CONST_G
 - ANGLE_DEGREE
 - ANGLE_ARCMIN
 - ANGLE_ARCSEC
 - ANGLE_MILLIARCSEC
 - ANGLE_MICROARCSEC

"""

import math

MASS_SUN          = 1.98855e30
"""Solar mass expressed in units of 1 kg."""

DIST_PC           = 3.0856775714409184e16
"""Distance of one parsec expressed in units of 1 m."""

DIST_KPC          = 1000.0*DIST_PC
"""Distance of one kiloparsec expressed in units of 1 m."""

DIST_MPC          = 1000.0*DIST_KPC
"""Distance of one megaparsec expressed in units of 1 m."""

DIST_GPC          = 1000.0*DIST_MPC
"""Distance of one gigaparsec expressed in units of 1 m."""

SPEED_C           = 299792458.0
"""The speed of light, expressed in units of 1 m/s."""

CONST_G           = 6.67300e-11
"""The gravitational constant, expressed in units of 1 m^3/(kg s)."""

ANGLE_DEGREE      = math.pi/180.0
"""Angle of one degree, expressed in radians."""

ANGLE_ARCMIN      = ANGLE_DEGREE/60.0
"""Angle one arcminute, expressed in radians."""

ANGLE_ARCSEC      = ANGLE_ARCMIN/60.0
"""Angle one arcsecond, expressed in radians."""

ANGLE_MILLIARCSEC = ANGLE_ARCSEC/1000.0
"""Angle one milliarcsecond, expressed in radians."""

ANGLE_MICROARCSEC = ANGLE_MILLIARCSEC/1000.0
"""Angle one microarcsecond, expressed in radians."""

