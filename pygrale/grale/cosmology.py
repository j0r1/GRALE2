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
            prod = abs(Wk)*s
            if abs(prod) < 1e-5: # flat enough
                A = s
            elif Wk > 0:
                A = 1.0/Wk**0.5 * math.sinh(prod)
            else: # Wk < 0
                A = 1.0/(-Wk)**0.5 * math.sin(prod)

        return (constants.SPEED_C*A*((1.0/(1.0+z2))/100000.0)/self.h)*constants.DIST_MPC

    def findRedshiftForAngularDiameterDistance(self, Dtarget, zref = 0, zmax = 20, iterations = 10, stepsPerIteration = 100):
        """This function attempts to find the redshift(s) that correspond to the
        angular diameter distance ``Dtarget``. 
        
        Parameters are:

         - ``Dtarget``: the angular diameter distance for which the possible redshifts
           should be searched.
         - ``zref``: in case the angular diameter distance should be calculated with
           respect to a specific redshift, this redshift can be specified here.
         - ``zmax``: the algorithm will look for redshifts that are smaller than this
           value.
         - ``iterations``: the procedure subdivides the possible range of redshift values
           into a number of pieces, yielding a first set of possible redshifts. These are
           then refined a number of times, specified by this value.
         - ``stepsPerIteration``: the number of parts a trial redshift range should be
           subdivided in.
        """

        zref, zmax = float(zref), float(zmax)

        def findIntervals(z1, z2, steps, iterations):

            intervals = [ ]
            dz = (z2-z1)/(steps-1)
            z = z1

            zD = [ ]

            # Calculate the mappings in this range
            for s in range(steps):
                D = self.getAngularDiameterDistance(zref, z)
                zD.append((z,D))
                z += dz

            # And see which intervals contain Dtarget
            for s in range(steps-1):
                if (zD[s][1] <= Dtarget and zD[s+1][1] > Dtarget) or (zD[s][1] >= Dtarget and zD[s+1][1] < Dtarget):
                    intervals.append((zD[s][0], zD[s+1][0]))

            #print(intervals)

            if iterations == 0:
                return intervals

            finerIntervals = [ ]

            for zf1, zf2 in intervals:
                sub = findIntervals(zf1, zf2, steps, iterations-1)
                if len(sub) >= 1:
                    finerIntervals += sub
                else:
                    finerIntervals += [( zf1, zf2)]

            return finerIntervals

        if zref < 0 or zmax < 0:
            raise CosmologyException("Both zref and zmax must be positive")
        if zref > zmax:
            raise CosmologyException("The value of zref must be less than zmax")

        intervals = findIntervals(zref, zmax, stepsPerIteration, iterations)
        zvalues = [ ]
        for z1, z2 in intervals:
            zvalues.append((z1+z2)/2.0)

        zvalues.sort()
        reducedZValues = []

        prevZ = -1e200
        for z in zvalues:
            if abs(z-prevZ) > 1e-6:
                reducedZValues.append(z)
            prevZ = z

        return reducedZValues

_defaultCosmology = [ None ]

def getDefaultCosmology():
    """TODO"""
    return _defaultCosmology[0]

def setDefaultCosmology(x):
    """TODO"""
    _defaultCosmology[0] = x

