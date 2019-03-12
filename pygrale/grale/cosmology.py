"""With this module you can calculate angular diameter distances for a specific
cosmology. Inversely, if you know the angular diameter distance to an object, you
can try to reconstruct the redshifts that are compatible with it.
"""
from __future__ import print_function
from . import constants
import scipy.integrate as igr
import math

class CosmologyException(Exception):
    """An exception that will be raised in case something goes wrong
    in the cosmology calculations."""
    pass

class Cosmology(object):
    """Specifies the cosmology parameters for which angular diameter distances
    will be calculated."""

    def __init__(self, h, Omega_m, Omega_r, Omega_v, w = -1.0):
        """
        Parameters are:

         - h: specifies the hubble constant :math:`H_0` as :math:`h \\times 100 \\; \\rm{km}\\; \\rm{s}^{-1} \\rm{Mpc}^{-1}`
         - Omega_m: matter density parameter :math:`\\Omega_m`
         - Omega_r: radiation density parameter :math:`\\Omega_r`
         - Omega_v: vacuum density parameter :math:`\\Omega_v`
         - w: specifies the equation of state of the vacuum energy as :math:`p = w \\times \\rho`

        """
        self.h = float(h)
        self.Wm = float(Omega_m)
        self.Wr = float(Omega_r)
        self.Wv = float(Omega_v)
        self.w = float(w)

    def getParameters(self):
        """Returns a dictionary that contains the model parameters that were specified in the
        constructor."""
        return { "h": self.h, "Omega_m": self.Wm, "Omega_r": self.Wr, "Omega_v": self.Wv, "w": self.w }

    def getAngularDiameterDistance(self, z1, z2 = None):
        """
        Returns the angular diameter distance between two redshifts, or if only one
        redshift is specified, between redshift 0 and the specified one. The returned distance
        is specified in meters, but can be converted using the constants in :mod:`grale.constants`.
        """
        if z1 is None:
            if z2 is None:
                raise CosmologyException("No redshifts given")
            else:
                z1, z2 = z2, z1
                
        if z2 is None:
            z1, z2 = 0.0, z1

        if z1 > z2:
            z1, z2 = z2, z1

        if z1 < 0 or z2 < 0:
            raise CosmologyException("Both redshifts must be positive")

        w = self.w
        Wm, Wr, Wv = self.Wm, self.Wr, self.Wv
        Wk = 1.0 - Wv - Wr - Wm

        def f(R):
            return 1.0/(Wm*R + Wr + Wv*R**(1.0-3.0*w) + Wk*R*R)**0.5
                    
        s = igr.quad(f, 1.0/(1.0+z2), 1.0/(1.0+z1))[0]
        if Wk == 0:
            A = s
        else:
            prod = abs(Wk)**0.5*s
            if abs(prod) < 1e-5: # flat enough
                A = s
            elif Wk > 0:
                A = 1.0/Wk**0.5 * math.sinh(prod)
            else: # Wk < 0
                A = 1.0/(-Wk)**0.5 * math.sin(prod)

        return (constants.SPEED_C*A*((1.0/(1.0+z2))/100000.0)/self.h)*constants.DIST_MPC

    def findRedshiftForAngularDiameterDistance(self, Dtarget, zref = 0, zmax = 20):
        """This function attempts to find the redshift(s) that correspond to the
        angular diameter distance ``Dtarget``. 
        
        Parameters are:

         - ``Dtarget``: the angular diameter distance for which the possible redshifts
           should be searched.
         - ``zref``: in case the angular diameter distance should be calculated with
           respect to a specific redshift, this redshift can be specified here.
         - ``zmax``: the algorithm will look for redshifts that are smaller than this
           value.
        """

        zref, zmax = float(zref), float(zmax)

        if zref < 0 or zmax < 0:
            raise CosmologyException("Both zref and zmax must be positive")
        if zref > zmax:
            raise CosmologyException("The value of zref must be less than zmax")

        import scipy.optimize as opt

        def f(z):
            if z < zref:
                z = zref
            if z > zmax:
                z = zmax
            return self.getAngularDiameterDistance(zref, z) - Dtarget

        z1 = opt.root(f, (zmax-zref)*0.01+zref)
        z2 = opt.root(f, (zmax-zref)*0.99+zref)
        
        sols = [ ]
        for z in z1, z2:
            if z.success: sols.append(float(z.x))
        return sols

_defaultCosmology = [ None ]

def getDefaultCosmology():
    """Returns the default cosmological model that has been set 
    using :func:`setDefaultCosmology`."""
    return _defaultCosmology[0]

def setDefaultCosmology(x):
    """Sets the default cosmological model to `x`, which should
    be an instance of :class:`Cosmology`."""
    _defaultCosmology[0] = x

