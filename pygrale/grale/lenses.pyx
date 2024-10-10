"""This module contains classes to work with different kinds of gravitational 
lens models.

The core class is :class:`GravitationalLens`, but this should not be instantiated
directly. Instead allocate a class derived from this; currently available are

 * :class:`GaussLens`
 * :class:`MultiplePlummerLens`
 * :class:`PlummerLens`
 * :class:`PointmassLens`
 * :class:`SISLens`
 * :class:`NSIELens`
 * :class:`NSISLens`
 * :class:`SIELens`
 * :class:`SquareLens`
 * :class:`MultipleSquareLens`
 * :class:`MultipleGaussLens`
 * :class:`MassSheetLens`
 * :class:`CompositeLens`
 * :class:`MassDiskLens`
 * :class:`PolynomialMassProfileLens`
 * :class:`TimeDelayAdjustLens`
 * :class:`ZeroMassLens`
 * :class:`MassDiskLensSmoothed`
 * :class:`MultipleWendlandLens`
 * :class:`DeflectionGridLens`
 * :class:`NFWLens`
 * :class:`EllipticNFWLens`
 * :class:`SersicLens`
 * :class:`EllipticSersicLens`
 * :class:`ProfileLens`
 * :class:`PIEMDLens`
 * :class:`PIMDLens`
 * :class:`AlphaPotLens`
 * :class:`HarmonicLens`
 * :class:`PotentialGridLens`
 * :class:`CircularPiecesLens`
 * :class:`MultiPlaneContainer`

"""

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cmap
from libcpp.memory cimport shared_ptr, unique_ptr, make_unique
from libcpp cimport bool
from libcpp.cast cimport dynamic_cast
from cython.operator cimport dereference as deref, postincrement
import cython
import numpy as np
cimport numpy as np
from cpython.array cimport array,clone
import struct
import random
from . import constants
from . import privutil
import copy
import multiprocessing

cimport grale.gravitationallens as gravitationallens
cimport grale.vector2d as vector2d
cimport grale.serut as serut
cimport grale.errut as errut
cimport grale.cppthreadslenscalc as threadcalc
cimport grale.imageplane as imageplane

include "stringwrappers.pyx"

ctypedef np.ndarray ndarray
ctypedef vector2d.Vector2Dd Vector2Dd

cdef int _numCalculationThreads = 0

experimentalThreads = False

class LensException(Exception):
    """This exception is raised if something goes wrong in the :class:`GravitationalLens` derived classes."""
    pass

def setNumberOfCalculationThreads(int n):
    """TODO"""
    if n < 0:
        raise LensException("Number of calculation threads must be at least zero (zero means all possible threads)")
    global _numCalculationThreads
    _numCalculationThreads = n

def getNumberOfCalculationThreads():
    """TODO"""
    global _numCalculationThreads
    if _numCalculationThreads == 0:
        _numCalculationThreads = multiprocessing.cpu_count()
    return _numCalculationThreads

# We use a random number as an identifier to make sure that GravitationalLens
# isn't allocated on its own, but only through subclassing
_gravLensRndId = int(random.random()*0xFFFFFFFF)

def getCriticalDensity(double Dd, double Ds, double Dds):
    """getCriticalDensity(Dd, Ds, Dds)
    
    Returns the critical density corresponding to angular diameter distances
    Dd (observer to lens), Ds (observer to source) and Dds (lens to source).

    This is defined as:

    .. code-block:: python

        Sigma_crit = c^2/(4*pi*G*Dd) * (Ds/Dds)

    """
    return constants.SPEED_C**2/(4.0*np.pi*constants.CONST_G*Dd)*(Ds/Dds)

cdef class GravitationalLens:
    """A base class for a gravitational lens."""

    cdef unique_ptr[gravitationallens.GravitationalLens] m_pLens

    cdef gravitationallens.GravitationalLens * _lens(self):
        return self.m_pLens.get()

    @staticmethod
    cdef gravitationallens.GravitationalLens * _getLens(GravitationalLens l):
        return l.m_pLens.get()

    @staticmethod
    cdef void _swapLenses(GravitationalLens l1, GravitationalLens l2):
        l1.m_pLens.swap(l2.m_pLens)

    def __init__(self, checkId = None):
        if checkId is None or checkId != _gravLensRndId:
            raise LensException("A 'GravitationalLens' object should not be allocated directly, only through subclasses.")

    cdef _check(self):
        if self._lens() == NULL:
            raise LensException("No internal lens object has been set")
    
    cdef _traceTheta1(self, double Ds, double Dds, np.ndarray[double, ndim=1] thetas):
        
        cdef Vector2Dd beta
        cdef np.ndarray[double,ndim=1] betas
        cdef int l,i

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            betas = np.zeros([l], dtype = np.double)
            for i in range(0,l,2):
                if not self._lens().traceTheta(Ds, Dds, Vector2Dd(thetas[i], thetas[i+1]), cython.address(beta)):
                    raise LensException(S(self._lens().getErrorString()))
                betas[i] = beta.getX()
                betas[i+1] = beta.getY()
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return betas

    cdef _getAlphaVector1(self, np.ndarray[double, ndim=1] thetas):
        
        cdef Vector2Dd alpha
        cdef np.ndarray[double,ndim=1] alphas
        cdef int l,i

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            alphas = np.zeros([l], dtype = np.double)
            for i in range(0,l,2):
                if not self._lens().getAlphaVector(Vector2Dd(thetas[i], thetas[i+1]), cython.address(alpha)):
                    raise LensException(S(self._lens().getErrorString()))
                alphas[i] = alpha.getX()
                alphas[i+1] = alpha.getY()
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return alphas

    cdef _getAlphaVectorDerivatives1(self, np.ndarray[double, ndim=1] thetas):
        
        cdef double axx = 0, ayy = 0, axy = 0
        cdef np.ndarray[double,ndim=1] derivatives
        cdef int l,i,j

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            derivatives = np.zeros([(l//2)*3], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                if not self._lens().getAlphaVectorDerivatives(Vector2Dd(thetas[i], thetas[i+1]), axx, ayy, axy):
                    raise LensException(S(self._lens().getErrorString()))
                derivatives[j] = axx
                derivatives[j+1] = ayy
                derivatives[j+2] = axy
                j += 3
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return derivatives

    cdef _getAlphaVectorSecondDerivatives1(self, np.ndarray[double, ndim=1] thetas):

        cdef double axxx = 0, ayyy = 0, axxy = 0, ayyx = 0
        cdef np.ndarray[double,ndim=1] derivatives
        cdef int l,i,j

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            derivatives = np.zeros([(l//2)*4], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                if not self._lens().getAlphaVectorSecondDerivatives(Vector2Dd(thetas[i], thetas[i+1]), axxx, ayyy, axxy, ayyx):
                    raise LensException(S(self._lens().getErrorString()))
                derivatives[j] = axxx
                derivatives[j+1] = ayyy
                derivatives[j+2] = axxy
                derivatives[j+3] = ayyx
                j += 4
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return derivatives

    cdef _getInverseMagnification1(self, double Ds, double Dds, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] invMags
        cdef int l,i,j

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            invMags = np.zeros([(l//2)], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                invMags[j] = self._lens().getInverseMagnification(Ds, Dds, Vector2Dd(thetas[i], thetas[i+1]))
                j += 1
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return invMags

    cdef _getSurfaceMassDensity1(self, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] dens
        cdef int l,i,j

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            dens = np.zeros([(l//2)], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                dens[j] = self._lens().getSurfaceMassDensity(Vector2Dd(thetas[i], thetas[i+1]))
                j += 1
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return dens

    cdef _getProjectedPotential1(self, double Ds, double Dds, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] potential
        cdef int l,i,j
        cdef double value = 0

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            potential = np.zeros([(l//2)], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                if not self._lens().getProjectedPotential(Ds, Dds, Vector2Dd(thetas[i], thetas[i+1]), cython.address(value)):
                    raise LensException(S(self._lens().getErrorString()))

                potential[j] = value
                j += 1
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return potential

    cdef _traceTheta1_new(self, double Ds, double Dds, np.ndarray[double, ndim=1] thetas):
        
        cdef Vector2Dd beta
        cdef np.ndarray[double,ndim=1] betas
        cdef int l,i
        cdef double *pThetaX
        cdef double *pThetaY
        cdef double *pBetaX
        cdef double *pBetaY
        cdef size_t thetaStride, betaStride
        cdef string errStr
        cdef size_t numTreads = getNumberOfCalculationThreads()

        print("_traceTheta1_new")

        if thetas.shape[0] % 2 == 0:
            pThetaX = &thetas[0]
            pThetaY = &thetas[1]
            
            l = thetas.shape[0]
            betas = np.zeros([l], dtype = np.double)
            pBetaX = &betas[0]
            pBetaY = &betas[1]

            thetaStride = 2
            betaStride = 2

            if not threadcalc.threadsTraceTheta(deref(self._lens()), errStr,
                    Ds, Dds, pThetaX, pThetaY, thetaStride, pBetaX, pBetaY, betaStride,
                    l//2, numTreads):
                raise LensException(S(errStr))

        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return betas

    cdef _getAlphaVector1_new(self, np.ndarray[double, ndim=1] thetas):
        
        cdef Vector2Dd alpha
        cdef np.ndarray[double,ndim=1] alphas
        cdef int l,i
        cdef double *pThetaX
        cdef double *pThetaY
        cdef double *pAlphaX
        cdef double *pAlphaY
        cdef size_t thetaStride, alphaStride
        cdef string errStr
        cdef size_t numTreads = getNumberOfCalculationThreads()

        print("_getAlphaVector1_new")

        if thetas.shape[0] % 2 == 0:
            pThetaX = &thetas[0]
            pThetaY = &thetas[1]

            l = thetas.shape[0]
            alphas = np.zeros([l], dtype = np.double)
            pAlphaX = &alphas[0]
            pAlphaY = &alphas[1]

            thetaStride = 2
            alphaStride = 2

            if not threadcalc.threadsGetAlphaVector(deref(self._lens()), errStr,
                    pThetaX, pThetaY, thetaStride, pAlphaX, pAlphaY, alphaStride,
                    l//2, numTreads):
                raise LensException(S(errStr))
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return alphas

    cdef _getAlphaVectorDerivatives1_new(self, np.ndarray[double, ndim=1] thetas):
        
        cdef double axx = 0, ayy = 0, axy = 0
        cdef np.ndarray[double,ndim=1] derivatives
        cdef int l,i,j
        cdef double *pThetaX
        cdef double *pThetaY
        cdef double *pAlphaXX
        cdef double *pAlphaYY
        cdef double *pAlphaXY
        cdef size_t thetaStride, alphaStride
        cdef string errStr
        cdef size_t numTreads = getNumberOfCalculationThreads()

        print("_getAlphaVectorDerivatives1_new")

        if thetas.shape[0] % 2 == 0:
            pThetaX = &thetas[0]
            pThetaY = &thetas[1]

            l = thetas.shape[0]
            derivatives = np.zeros([(l//2)*3], dtype = np.double)
            pAlphaXX = &derivatives[0]
            pAlphaYY = &derivatives[1]
            pAlphaXY = &derivatives[2]

            thetaStride = 2
            alphaStride = 3

            if not threadcalc.threadsGetAlphaVectorDerivatives(deref(self._lens()), errStr,
                    pThetaX, pThetaY, thetaStride, pAlphaXX, pAlphaYY, pAlphaXY, alphaStride,
                    l//2, numTreads):
                raise LensException(S(errStr))
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return derivatives

    cdef _getAlphaVectorSecondDerivatives1_new(self, np.ndarray[double, ndim=1] thetas):

        cdef double axxx = 0, ayyy = 0, axxy = 0, ayyx = 0
        cdef np.ndarray[double,ndim=1] derivatives
        cdef int l,i,j
        cdef double *pThetaX
        cdef double *pThetaY
        cdef double *pAlphaXXX
        cdef double *pAlphaYYY
        cdef double *pAlphaXXY
        cdef double *pAlphaYYX
        cdef size_t thetaStride, alphaStride
        cdef string errStr
        cdef size_t numTreads = getNumberOfCalculationThreads()

        print("_getAlphaVectorSecondDerivatives1_new")

        if thetas.shape[0] % 2 == 0:
            pThetaX = &thetas[0]
            pThetaY = &thetas[1]

            l = thetas.shape[0]
            derivatives = np.zeros([(l//2)*4], dtype = np.double)
            pAlphaXXX = &derivatives[0]
            pAlphaYYY = &derivatives[1]
            pAlphaXXY = &derivatives[2]
            pAlphaYYX = &derivatives[3]

            thetaStride = 2
            alphaStride = 4

            if not threadcalc.threadsGetAlphaVectorSecondDerivatives(deref(self._lens()), errStr,
                    pThetaX, pThetaY, thetaStride, pAlphaXXX, pAlphaYYY, pAlphaXXY, pAlphaYYX,
                    alphaStride, l//2, numTreads):
                raise LensException(S(errStr))
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return derivatives


    cdef _getInverseMagnification1_new(self, double Ds, double Dds, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] invMags
        cdef int l,i,j
        cdef double *pThetaX
        cdef double *pThetaY
        cdef double *pInvMag
        cdef size_t thetaStride, invStride
        cdef string errStr
        cdef size_t numTreads = getNumberOfCalculationThreads()

        print("_getInverseMagnification1_new")

        if thetas.shape[0] % 2 == 0:
            pThetaX = &thetas[0]
            pThetaY = &thetas[1]
            
            l = thetas.shape[0]
            invMags = np.zeros([(l//2)], dtype = np.double)
            pInvMag = &invMags[0]

            thetaStride = 2
            invStride = 1
            
            if not threadcalc.threadsGetInverseMagnification(deref(self._lens()), errStr,
                    Ds, Dds, pThetaX, pThetaY, thetaStride, pInvMag, invStride,
                    l//2, numTreads):
                raise LensException(S(errStr))
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return invMags

    cdef _getSurfaceMassDensity1_new(self, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] dens
        cdef int l,i,j
        cdef double *pThetaX
        cdef double *pThetaY
        cdef double *pDens
        cdef size_t thetaStride, densStride
        cdef string errStr
        cdef size_t numTreads = getNumberOfCalculationThreads()

        print("_getSurfaceMassDensity1_new")

        if thetas.shape[0] % 2 == 0:
            pThetaX = &thetas[0]
            pThetaY = &thetas[1]

            l = thetas.shape[0]
            dens = np.zeros([(l//2)], dtype = np.double)
            
            pDens = &dens[0]

            thetaStride = 2
            densStride = 1
            
            if not threadcalc.threadsGetSurfaceMassDensity(deref(self._lens()), errStr,
                    pThetaX, pThetaY, thetaStride, pDens, densStride,
                    l//2, numTreads):
                raise LensException(S(errStr))
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return dens

    cdef _getProjectedPotential1_new(self, double Ds, double Dds, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] potential
        cdef int l,i,j
        cdef double value = 0
        cdef double *pThetaX
        cdef double *pThetaY
        cdef double *pPot
        cdef size_t thetaStride, potStride
        cdef string errStr
        cdef size_t numTreads = getNumberOfCalculationThreads()

        print("_getProjectedPotential1_new")

        if thetas.shape[0] % 2 == 0:
            pThetaX = &thetas[0]
            pThetaY = &thetas[1]
            
            l = thetas.shape[0]
            potential = np.zeros([(l//2)], dtype = np.double)
            
            pPot = &potential[0]

            thetaStride = 2
            potStride = 1

            if not threadcalc.threadsGetProjectedPotential(deref(self._lens()), errStr,
                    Ds, Dds, pThetaX, pThetaY, thetaStride, pPot, potStride,
                    l//2, numTreads):
                raise LensException(S(errStr))
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return potential

    cdef _getRadialMassProfile1(self, np.ndarray[double, ndim=1] thetaRadii):

        cdef int i,l
        cdef double value = 0
        cdef np.ndarray[double,ndim=1] profile
        cdef gravitationallens.SymmetricLens *pSymLens = gravitationallens.SymmetricLens.cast(self._lens())
        if not pSymLens:
            raise LensException("The lens needs to be a circularly symmetric lens type for this to work")

        l = thetaRadii.shape[0]
        profile = np.zeros([l], dtype = np.double)
        for i in range(l):
            profile[i] = pSymLens.getMassInside(thetaRadii[i])

        return profile

    cdef _getRadialDensityProfile1(self, np.ndarray[double, ndim=1] thetaRadii):

        cdef int i,l
        cdef double value = 0
        cdef np.ndarray[double,ndim=1] profile
        cdef gravitationallens.SymmetricLens *pSymLens = gravitationallens.SymmetricLens.cast(self._lens())
        if not pSymLens:
            raise LensException("The lens needs to be a circularly symmetric lens type for this to work")

        l = thetaRadii.shape[0]
        profile = np.zeros([l], dtype = np.double)
        for i in range(l):
            profile[i] = pSymLens.getProfileSurfaceMassDensity(thetaRadii[i])

        return profile

    cdef _reshapeAndCall1D(self, functionName, thetas, int coreNumIn, int coreNumOut):
        cdef int totalElements, l, i

        self._check()

        l = len(thetas.shape)
        if l == 1:
            return functionName(thetas)

        if l > 1 and thetas.shape[l-1] == coreNumIn:
            totalElements = 1
            for i in range(l):
                totalElements *= thetas.shape[i]

            outShape = thetas.shape[:-1] 
            if coreNumOut > 1:
                outShape += (coreNumOut,)
            return np.reshape(functionName(np.reshape(thetas,[totalElements])), outShape)

        raise LensException("Bad array dimensions")

    def traceTheta(self, Ds, Dds, thetas):
        """traceTheta(Ds, Dds, thetas)
        
        Traces the image plane positions in 'thetas' to the source plane that corresponds
        to the angular diameter distances Ds and Dds.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._traceTheta1_new(float(Ds), float(Dds), x) if experimentalThreads else self._traceTheta1(float(Ds), float(Dds), x), thetas, 2, 2)

    def getAlphaVector(self, thetas):
        """getAlphaVector(thetas)
        
        Returns the deflection angles at the positions in 'thetas'.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getAlphaVector1_new(x) if experimentalThreads else self._getAlphaVector1(x), thetas, 2, 2)

    def getAlphaVectorDerivatives(self, thetas):
        r"""getAlphaVectorDerivatives(thetas)
        
        Returns the derivatives of the deflection angle at the positions specified in 'thetas'.
        For each position (thetax, thetay), three values (axx, ayy, axy) are returned, where

        .. math::

            \begin{array}{l}
              {\rm axx} = \frac{\partial \hat{\alpha}_x}{\partial \theta_x} \\
              {\rm ayy} = \frac{\partial \hat{\alpha}_y}{\partial \theta_y} \\
              {\rm axy} = \frac{\partial \hat{\alpha}_x}{\partial \theta_y} = \frac{\partial \hat{\alpha}_y}{\partial \theta_x} 
            \end{array}

        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getAlphaVectorDerivatives1_new(x) if experimentalThreads else self._getAlphaVectorDerivatives1(x), thetas, 2, 3)

    def getAlphaVectorSecondDerivatives(self, thetas):
        r"""getAlphaVectorSecondDerivatives(thetas)

        Returns the second derivatives of the deflection angle at the positions specified in 'thetas'.
        For each position (thetax, thetay), three values (axxx, ayyy, axxy, ayyx) are returned, where

        TODO
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getAlphaVectorSecondDerivatives1_new(x) if experimentalThreads else self._getAlphaVectorSecondDerivatives1(x), thetas, 2, 4)

    def getInverseMagnification(self, Ds, Dds, thetas):
        """getInverseMagnification(Ds, Dds, thetas)
        
        Returns the inverse magnification values at the positions in 'thetas', for a source
        plane with angular diameter distances Ds and Dds.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getInverseMagnification1_new(float(Ds), float(Dds), x) if experimentalThreads else self._getInverseMagnification1(float(Ds), float(Dds), x), thetas, 2, 1)

    def getSurfaceMassDensity(self, thetas):
        """getSurfaceMassDensity(thetas)
        
        Returns the surface mass density of the gravitational lens at the positions
        in 'thetas'.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getSurfaceMassDensity1_new(x) if experimentalThreads else self._getSurfaceMassDensity1(x), thetas, 2, 1)

    def getSurfaceMassDensityMap(self, bottomLeft, topRight, int numX, int numY, feedbackObject = "default", renderer = "default", reduceToPixels = True):
        """getSurfaceMassDensityMap(bottomLeft, topRight, numX, numY, feedbackObject = "default", renderer = "default", reduceToPixels = True)
        
        Returns the surface density map for the area between the bottomLeft and topRight
        coordinates. The surface mass density is sampled at numX and numY points and
        these values are returned directly if reduceToPixels is False. The grid can
        also be seen as a grid of pixels, and if reduceToPixels is True, the values will
        be averaged to obtain (numX-1)*(numY-1) pixel values.

        The feedbackObject specifies the method which is used for feedback when
        calculating the values as the grid points. See :mod:`grale.feedback` for
        more information.

        If renderer is None (the default), calculations will occur in the same process,
        on a single CPU core. For other possible values, see the :mod:`grale.renderers`
        module.
        """
        cdef double width = topRight[0] - bottomLeft[0]
        cdef double height = topRight[1] - bottomLeft[1]
        cdef double x0, y0, dX, dY
        cdef np.ndarray[double,ndim=2] massMap
        cdef int xi, yi, pos, advance
        cdef double value, pct

        self._check()

        if width <= 0 or height <= 0:
            raise LensException("Both width and height of the region must be positive")
        if numX < 2 or numY < 2:
            raise LensException("Number of points in X and Y direction must be at least 2")

        # Makes sure the feedback object of the renderer matches the one specified (if any), and returns
        # the used one, or a dummy
        renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "MASSDENS")

        massMap = np.zeros([numY, numX], dtype = np.double)

        if renderer is None:
            x0 = bottomLeft[0]
            y0 = bottomLeft[1]
            dX = width/(numX-1)
            dY = height/(numY-1)

            feedbackObject.onStatus("Calculating mass map")

            for yi in range(numY):
                y = y0 + dY*yi

                pct = 100.0 * <double>(yi)/<double>(numY) + 0.5
                feedbackObject.onProgress(<int>(pct))

                for xi in range(numX):
                    x = x0 + dX*xi
                    massMap[yi,xi] = self._lens().getSurfaceMassDensity(Vector2Dd(x,y))

            feedbackObject.onProgress(100)
            feedbackObject.onStatus("Done")

        else:
            if renderer.renderType != "MASSDENS":
                raise LensException("Specified renderer does not appear to be for the MASSDENS type")

            massBytes = renderer.render(self.toBytes(), bottomLeft, topRight, numX, numY)
            pos = 0
            advance = struct.calcsize("=d")

            for yi in range(numY):
                for xi in range(numX):
                    value = struct.unpack_from("=d", massBytes, pos)[0]
                    pos += advance
                    massMap[yi,xi] = value

        if reduceToPixels:
            massMap = 0.25*massMap[0:numY-1,0:numX-1] + 0.25*massMap[0:numY-1,1:numX] + 0.25*massMap[1:numY,0:numX-1] + 0.25*massMap[1:numY,1:numX]

        return massMap

    def getProjectedPotential(self, Ds, Dds, thetas):
        """getProjectedPotential(Ds, Dds, thetas)
        
        Returns the projected potential for a source plane with angular diameter
        distances Ds and Dds, measured at the positions in 'thetas'.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getProjectedPotential1_new(float(Ds), float(Dds), x) if experimentalThreads else self._getProjectedPotential1(float(Ds), float(Dds), x), thetas, 2, 1)

    def getRadialMassProfile(self, thetaRadii):
        """getRadialMassProfile(thetaRadii)

        If the lens is a symmetric lens, this function returns the radial mass
        profile at the specified radii.
        """
        thetaRadii = np.array(thetaRadii)
        return self._reshapeAndCall1D(lambda x : self._getRadialMassProfile1(x), thetaRadii, 1, 1)

    def getRadialDensityProfile(self, thetaRadii):
        """getRadialDensityProfile(thetaRadii)

        If the lens is a symmetric lens, this function returns the density
        profile at the specified radii.
        """
        thetaRadii = np.array(thetaRadii)
        return self._reshapeAndCall1D(lambda x : self._getRadialDensityProfile1(x), thetaRadii, 1, 1)

    def getLensDistance(self):
        """getLensDistance()
        
        Returns the angular diameter distance Dd that was specified when creating
        this gravitational lens object.
        """
        self._check()
        return self._lens().getLensDistance()

    def setLensDistance(self, double Dd):
        """setLensDistance(Dd)
        
        This function changes the angular diameter distance to the lens to Dd.
        """
        self._check()
        if Dd <= 0:
            raise LensException("Invalid distance")
        self._lens().setLensDistance(Dd)

    def getTimeDelay(self, zd, Ds, Dds, theta, beta):
        r"""getTimeDelay(zd, Ds, Dds, thetas, beta)
        
        Returns the time delay (up to an unknown constant) for this gravitational
        lens with redshift 'zd', for a source position 'beta' in a source plane
        specified by the angular diameter distances 'Ds' and 'Dds', and for an
        image plane positions 'thetas'.

        Note that for this function, the 'thetas' positions do not necessarily
        need to trace back to 'beta'. This allows, for example, to calculate the
        time delays for a grid of image plane points (and one fixed 'beta') and 
        identify the locations where images are seen as the points where the
        time delay reaches a local optimum.
        """
        zd, Ds, Dds = float(zd), float(Ds), float(Dds)
        if zd < 0 or Ds <= 0 or Dds <= 0:
            raise LensException("Bad distance or redshift")

        potential = self.getProjectedPotential(Ds, Dds, theta)
        diff = theta - beta
        
        totalElements = 1
        for i in theta.shape:
            totalElements *= i

        diff2 = diff*diff
        diff2 = np.reshape(diff2, [totalElements//2, 2])
        distances = np.reshape(diff2[:,0] + diff2[:,1], potential.shape)
        factor = self._lens().getLensDistance() * (1.0 + zd) / constants.SPEED_C * (Ds/Dds)
        delay = (0.5*distances - potential)*factor
        return delay

    def getSuggestedScales(self):
        """getSuggestedScales()

        Returns a dictionary with suggested deflection and potential scales, which
        can be used in the OpenCL routines to be able to calculate in single precision
        floating point."""
        cdef double defScale = 0, potScale = 0

        self._check()
        if not self._lens().getSuggestedScales(&defScale, &potScale):
            raise LensException(S(self._lens().getErrorString()))
        return { "deflectionscale": defScale, "potentialscale": potScale }

    def getCLParameterCounts(self):
        """getCLParameterCounts()

        Returns a tuple containing the number of integer parameters and floating 
        point parameters this lens's OpenCL implementation needs."""
        cdef int numIntParams = 0, numFloatParams = 0

        self._check()
        if not self._lens().getCLParameterCounts(&numIntParams, &numFloatParams):
            raise LensException(S(self._lens().getErrorString()))
        return (numIntParams, numFloatParams)

    def getCLParameters(self, deflectionscale, potentialscale):
        """getCLParameters(deflectionscale, potentialscale)

        Returns a tuple of integer and floating point parameters to use this lens in OpenCL,
        expressed using the specified deflection and potential scales."""

        cdef np.ndarray[np.float32_t,ndim=1] floatParams
        cdef np.ndarray[np.int32_t,ndim=1] intParams

        numIntParams, numFloatParams = self.getCLParameterCounts()
        intParams = np.empty((numIntParams,),dtype=np.int32)
        floatParams = np.empty((numFloatParams,),dtype=np.float32)
        if not self._lens().getCLParameters(deflectionscale, potentialscale, <int*>intParams.data, <float*>floatParams.data):
            raise LensException(S(self._lens().getErrorString()))
        return (intParams, floatParams)

    def getCLLensProgram(self, deflectionScale, potentialScale, derivatives=True, potential=True):
        """getCLLensProgram(deflectionScale, potentialScale, derivatives=True, potential=True)

        Returns a concatenation of :func:`getCLLensQuantitiesStructure` and
        :func:`getCLProgram`"""
        cdef string subName
        
        self._check()
        program = self._lens().getCLLensProgram(deflectionScale, potentialScale, subName, derivatives, potential)
        return (S(subName), S(program))

    def getCLProgram(self, deflectionScale, potentialScale, derivatives=True, potential=True):
        """getCLProgram(deflectionScale, potentialScale, derivatives=True, potential=True)

        Returns an OpenCL program to calculate the deflection and optionally
        deflection derivatives and lensing potential, for this lens model."""
        cdef string subName

        self._check() 
        program = self._lens().getCLProgram(deflectionScale, potentialScale, subName, derivatives, potential)
        return (S(subName), S(program))

    def getCLLensQuantitiesStructure(self, derivatives=True, potential=True):
        """getCLLensQuantitiesStructure(derivatives=True, potential=True)

        The OpenCL routine returned by :func:`getCLProgram` returns a
        structure, for which the code can be queried by this function.
        """
        self._check()
        return S(self._lens().getCLLensQuantitiesStructure(derivatives, potential))

    def setDerivativeAngularDistanceScale(self, double scale):
        """setDerivativeAngularDistanceScale(scale)
        
        In case the :func:`getAlphaVectorDerivatives` function is not provided
        for a specific lens, the derivatives will be calculated numerically
        using points that are this distance apart.
        """
        self._check()
        self._lens().setDerivativeAngularDistanceScale(scale)

    def getCriticalDensity(self, double Ds, double Dds):
        """getCriticalDensity(Ds, Dds)
        
        Calls :func:`getCriticalDensity` where Dd is the angular diameter
        distance of this gravitational lens (see :func:`getLensDistance`)
        """
        Dd = self.getLensDistance()
        return getCriticalDensity(Dd, Ds, Dds)

    @staticmethod
    def load(fileName):
        """load(fileName)
        
        Loads the gravitational lens that was previously saved (see :func:`save`) to the file
        with name 'fileName'.
        """
        cdef string errorString
        cdef unique_ptr[gravitationallens.GravitationalLens] pLens
        if not gravitationallens.GravitationalLens.load(B(fileName), pLens, errorString):
            raise LensException(S(errorString))

        return GravitationalLens._finalizeLoadedLens(pLens)

    @staticmethod
    cdef _finalizeLoadedLens(unique_ptr[gravitationallens.GravitationalLens] &pLens):
        cdef gravitationallens.LensType t = deref(pLens).getLensType()
        if t == gravitationallens.Gaussian:
            l = GaussLens(None, None)
        elif t == gravitationallens.MultiplePlummers:
            l = MultiplePlummerLens(None, None)
        elif t == gravitationallens.Plummer:
            l = PlummerLens(None, None)
        elif t == gravitationallens.Pointmass:
            l = PointmassLens(None, None)
        elif t == gravitationallens.SIS:
            l = SISLens(None, None)
        elif t == gravitationallens.NSIE:
            l = NSIELens(None, None)
        elif t == gravitationallens.NSIS:
            l = NSISLens(None, None)
        elif t == gravitationallens.SIE:
            l = SIELens(None, None)
        elif t == gravitationallens.Square:
            l = SquareLens(None, None)
        elif t == gravitationallens.MultipleSquares:
            l = MultipleSquareLens(None, None)
        elif t == gravitationallens.MultipleGaussians:
            l = MultipleGaussLens(None, None)
        elif t == gravitationallens.MassSheet:
            l = MassSheetLens(None, None)
        elif t == gravitationallens.Composite:
            l = CompositeLens(None, None)
        elif t == gravitationallens.MassDisk:
            l = MassDiskLens(None, None)
        elif t == gravitationallens.Profile:
            l = ProfileLens(None, None)
        elif t == gravitationallens.PolynomialMassProfile:
            l = PolynomialMassProfileLens(None, None)
        elif t == gravitationallens.MultipleWendland:
            l = MultipleWendlandLens(None, None)
        elif t == gravitationallens.DeflectionGrid:
            l = DeflectionGridLens(None, None)
        elif t == gravitationallens.NFW:
            l = NFWLens(None, None)
        elif t == gravitationallens.EllipticNFW:
            l = EllipticNFWLens(None, None)
        elif t == gravitationallens.Sersic:
            l = SersicLens(None, None)
        elif t == gravitationallens.EllipticSersic:
            l = EllipticSersicLens(None, None)
        elif t == gravitationallens.PIEMD:
            l = PIEMDLens(None, None)
        elif t == gravitationallens.PIMD:
            l = PIMDLens(None, None)
        elif t == gravitationallens.AlphaPot:
            l = AlphaPotLens(None, None)
        elif t == gravitationallens.Harmonic:
            l = HarmonicLens(None, None)
        elif t == gravitationallens.PotentialGrid:
            l = PotentialGridLens(None, None)
        elif t == gravitationallens.CircularPieces:
            l = CircularPiecesLens(None, None)
        elif t == gravitationallens.MPContainer:
            l = MultiPlaneContainer(None, None)
        else: # Unknown, can still use the interface
            l = GravitationalLens(_gravLensRndId)

        l.m_pLens.swap(pLens)
        return l

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(bytes b)
        
        Instantiates a gravitational lens object from the byte array previously generated
        by :func:`toBytes`.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef unique_ptr[serut.MemorySerializer] m = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(b), <void*>NULL, 0)
        cdef unique_ptr[gravitationallens.GravitationalLens] pLens
        cdef string errorString
        
        if not gravitationallens.GravitationalLens.read(deref(m), pLens, errorString):
            raise LensException(S(errorString))

        return GravitationalLens._finalizeLoadedLens(pLens)

    def save(self, fileName):
        """save(fileName)
        
        Write this gravitational lens object to the file with name 'fileName'.
        Can be loaded again using :func:`load`.
        """
        self._check()
        if not self._lens().save(B(fileName)):
            raise LensException(S(self._lens().getErrorString()))

    def refinePosition(self, double Ds, double Dds, beta, startTheta, int numIterations = 4):
        cdef Vector2Dd vBeta = Vector2Dd(beta[0], beta[1])
        cdef Vector2Dd vStartTheta = Vector2Dd(startTheta[0], startTheta[1])
        cdef Vector2Dd theta = vStartTheta
        cdef string errStr

        self._check()
        if not imageplane.ImagePlane.staticRefinePosition(deref(self._lens()), Ds, Dds, vBeta, vStartTheta, theta, numIterations, errStr):
            raise LensException(S(errStr))
        return np.array([theta.getX(), theta.getY()], dtype=np.double)

    def toBytes(self):
        """toBytes()
        
        Converts the lens object to a byte array, which can be converted to a real
        lens object again using :func:`fromBytes`
        """
        cdef serut.VectorSerializer vSer

        self._check()
        if not self._lens().write(vSer):
            raise LensException(S(self._lens().getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

    # Will be overridden
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return NULL

    # Will be overridden
    cdef gravitationallens.GravitationalLensParams *_allocParams(self, params) except NULL:
        return NULL

    def _lensInit(self, Dd, params):
        cdef unique_ptr[gravitationallens.GravitationalLensParams] lp
        cdef unique_ptr[gravitationallens.GravitationalLens] pLens

        if Dd is None and params is None: # Don't initialize yet
            return

        lp = unique_ptr[gravitationallens.GravitationalLensParams](self._allocParams(params))
        pLens = unique_ptr[gravitationallens.GravitationalLens](self._allocLens())
        
        if not deref(pLens).init(Dd, lp.get()):
            err = deref(pLens).getErrorString()
            raise LensException(S(err))
        
        self.m_pLens.swap(pLens)

    @staticmethod
    def _checkParams(params, names):
        for n in names:
            if not n in params:
                raise LensException("Name '{}' is not present in the lens parameters".format(n))

        names = set(names)
        for n in params:
            if not n in names:
                raise LensException("Name '{}' in the lens parameters is not recognized".format(n))

    def getLensParameters(self):
        """getLensParameters()

        Obtain the parameters (or equivalent) that were used to create this gravitational
        lens instance.
        """
        raise LensException("Re-analyzing parameters has not been implemented for this lens type")

    def __getstate__(self):
        return self.toBytes()

    def __setstate__(self, state):
        l = GravitationalLens.fromBytes(state)
        GravitationalLens._swapLenses(self, l)

cdef class GaussLens(GravitationalLens):
    """A gravitational lens with a gaussian profile"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.GaussLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["mass", "width"])
        return new gravitationallens.GaussLensParams(params["mass"], params["width"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'mass':  the total mass of the lens
           * 'width': the (angular) width of the lens
        """
        super(GaussLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.GaussLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.GaussLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a GaussLens")

        return { "mass": pParams.getMass(), "width": pParams.getAngularWidth() }

cdef class MultiplePlummerLens(GravitationalLens):
    """A gravitational lenses consisting of a number of Plummer basis functions"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.MultiplePlummerLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef vector[gravitationallens.PlummerLensInfo] lensInfo

        for p in params:
            GravitationalLens._checkParams(p, ["mass", "width", "x", "y"])
            lensInfo.push_back(gravitationallens.PlummerLensInfo(p["mass"], p["width"], Vector2Dd(p["x"], p["y"])))

        return new gravitationallens.MultiplePlummerLensParams(lensInfo)

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a list of which each entry is a dictionary with the following entries:

           * 'mass':  the total mass of this Plummer basis function
           * 'width': the (angular) width of this basis function
           * 'x':     the X-position of the center of this basis function
           * 'y':     the Y-position of the center of this basis function
        """
        super(MultiplePlummerLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)
    
    def getLensParameters(self):
        cdef gravitationallens.MultiplePlummerLensParamsPtrConst pParams
        cdef vector[gravitationallens.PlummerLensInfo] lensInfo
        cdef int i, l
        
        self._check()
        pParams = dynamic_cast[gravitationallens.MultiplePlummerLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a MultiplePlummerLens")

        params = [ ]
        lensInfo = pParams.getLensInfo()
        l = lensInfo.size()
        for i in range(l):
            params.append({
                "mass": lensInfo[i].getMass(),
                "width": lensInfo[i].getAngularWidth(),
                "x": lensInfo[i].getAngularPosition().getX(),
                "y": lensInfo[i].getAngularPosition().getY()
            })

        return params

cdef class PlummerLens(GravitationalLens):
    """A gravitational lens corresponding to a projected Plummer sphere"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PlummerLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["mass", "width"])
        return new gravitationallens.PlummerLensParams(params["mass"], params["width"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'mass':  the total mass of this Plummer lens
           * 'width': the (angular) width of the Plummer profile
        """
        super(PlummerLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)
    
    def getLensParameters(self):
        cdef gravitationallens.PlummerLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.PlummerLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a PlummerLens")

        return { "mass": pParams.getLensMass(), "width": pParams.getAngularWidth() }

cdef class PointmassLens(GravitationalLens):
    """A point mass gravitational lens"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PointmassLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["mass"])
        return new gravitationallens.PointmassLensParams(params["mass"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'mass':  the total mass of this point mass lens
        """
        super(PointmassLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.PointmassLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.PointmassLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a PointmassLens")

        return { "mass": pParams.getLensMass() }

cdef class SISLens(GravitationalLens):
    """A singular isothermal sphere (SIS) gravitational lens"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.SISLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["velocityDispersion"])
        return new gravitationallens.SISLensParams(params["velocityDispersion"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'velocityDispersion': the velocity dispersion of the SIS lens
        """
        super(SISLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.SISLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.SISLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a SISLens")

        return { "velocityDispersion": pParams.getVelocityDispersion() }

cdef class NSIELens(GravitationalLens):
    """A non-singular isothermal ellipse (NSIE) gravitational lens"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.NSIELens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["velocityDispersion", "ellipticity", "coreRadius"])
        return new gravitationallens.NSIELensParams(params["velocityDispersion"], params["ellipticity"], params["coreRadius"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'velocityDispersion': the velocity dispersion of the NSIE lens
           * 'ellipticity':        the ellipticity parameter for this lens
           * 'coreRadius':         the (angular) radius of the core of the NSIE profile
        """
        super(NSIELens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.NSIELensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.NSIELensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a NSIELens")

        return { 
            "velocityDispersion": pParams.getVelocityDispersion(),
            "ellipticity": pParams.getEllipticity(),
            "coreRadius": pParams.getAngularCoreRadius()
        }

cdef class NSISLens(GravitationalLens):
    """A non-singular isothermal sphere (NSIS) gravitational lens"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.NSISLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["velocityDispersion", "coreRadius"])
        return new gravitationallens.NSISLensParams(params["velocityDispersion"], params["coreRadius"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'velocityDispersion': the velocity dispersion of the NSIS lens
           * 'coreRadius':         the (angular) radius of the core of the NSIS profile
        """
        super(NSISLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.NSISLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.NSISLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a NSISLens")

        return { 
            "velocityDispersion": pParams.getVelocityDispersion(),
            "coreRadius": pParams.getAngularCoreRadius()
        }

cdef class SIELens(GravitationalLens):
    r"""A singular isothermal ellipse (SIE) gravitational lens

.. math::

    \Sigma(\vec{\theta}) = \frac{\sigma_v^2\sqrt{f}}{2 G D_d\sqrt{\theta_x^2+f^2\theta_y^2}}

    \vec{\hat{\alpha}}(\vec{\theta}) = \frac{4 \pi \sigma_v^2}{c^2} \frac{\sqrt{f}}{\sqrt{1-f^2}}
 					   \left[
  	                                       \mathrm{asinh}\left(\frac{\sqrt{1-f^2}}{f}\frac{\theta_x}{|\vec{\theta}|}\right)\vec{e}_x +
  	                                       \mathrm{asin}\left(\sqrt{1-f^2}\frac{\theta_y}{|\vec{\theta}|}\right)\vec{e}_y
					   \right]

**References**

 *  `Kormann, R., Schneider, P., Bartelmann, M., Isothermal elliptical gravitational lens models. Astronomy and Astrophysics, 284:285-299, April 1994. <http://adsabs.harvard.edu/abs/1994A&A...284..285K>`_

"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.SIELens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["velocityDispersion", "ellipticity"])
        return new gravitationallens.SIELensParams(params["velocityDispersion"], params["ellipticity"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'velocityDispersion': the velocity dispersion of the SIE lens
           * 'ellipticity':        the ellipticity parameter for this lens
        """
        super(SIELens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.SIELensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.SIELensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a SIELens")

        return { 
            "velocityDispersion": pParams.getVelocityDispersion(),
            "ellipticity": pParams.getEllipticity(),
        }

cdef class SquareLens(GravitationalLens):
    """The gravitational lens effect of a square shaped region of contant surface mass density"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.SquareLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, ["mass", "width"])
        return new gravitationallens.SquareLensParams(params["mass"], params["width"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'mass':  the total mass of this square shaped lens
           * 'width': the (angular) width of the lens
        """
        super(SquareLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.SquareLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.SquareLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a SquareLens")

        return { "mass": pParams.getLensMass(), "width": pParams.getAngularWidth() }

cdef class MultipleSquareLens(GravitationalLens):
    """A gravitational lenses consisting of a number of square shaped basis functions"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.MultipleSquareLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef vector[gravitationallens.SquareLensInfo] lensInfo

        for p in params:
            GravitationalLens._checkParams(p, [ "mass", "width", "x", "y"])
            lensInfo.push_back(gravitationallens.SquareLensInfo(p["mass"], p["width"], Vector2Dd(p["x"], p["y"])))

        return new gravitationallens.MultipleSquareLensParams(lensInfo)

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a list of which each entry is a dictionary with the following entries:

           * 'mass':  the total mass of this square basis function
           * 'width': the (angular) width of this basis function
           * 'x':     the X-position of the center of this basis function
           * 'y':     the Y-position of the center of this basis function
        """
        super(MultipleSquareLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.MultipleSquareLensParamsPtrConst pParams
        cdef vector[gravitationallens.SquareLensInfo] lensInfo
        cdef int i, l
        
        self._check()
        pParams = dynamic_cast[gravitationallens.MultipleSquareLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a MultipleSquareLens")

        params = [ ]
        lensInfo = pParams.getLensInfo()
        l = lensInfo.size()
        for i in range(l):
            params.append({
                "mass": lensInfo[i].getMass(),
                "width": lensInfo[i].getAngularWidth(),
                "x": lensInfo[i].getAngularPosition().getX(),
                "y": lensInfo[i].getAngularPosition().getY()
            })

        return params

cdef class MultipleGaussLens(GravitationalLens):
    """A gravitational lenses consisting of a number of gaussian basis functions"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.MultipleGaussLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef vector[gravitationallens.GaussLensInfo] lensInfo

        for p in params:
            GravitationalLens._checkParams(p, ["mass", "width", "x", "y"])
            lensInfo.push_back(gravitationallens.GaussLensInfo(p["mass"], p["width"], Vector2Dd(p["x"], p["y"])))

        return new gravitationallens.MultipleGaussLensParams(lensInfo)

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a list of which each entry is a dictionary with the following entries:

           * 'mass':  the total mass of this gaussian basis function
           * 'width': the (angular) width of this basis function
           * 'x':     the X-position of the center of this basis function
           * 'y':     the Y-position of the center of this basis function
        """
        super(MultipleGaussLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.MultipleGaussLensParamsPtrConst pParams
        cdef vector[gravitationallens.GaussLensInfo] lensInfo
        cdef int i, l
        
        self._check()
        pParams = dynamic_cast[gravitationallens.MultipleGaussLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a MultipleGaussLens")

        params = [ ]
        lensInfo = pParams.getLensInfo()
        l = lensInfo.size()
        for i in range(l):
            params.append({
                "mass": lensInfo[i].getMass(),
                "width": lensInfo[i].getAngularWidth(),
                "x": lensInfo[i].getAngularPosition().getX(),
                "y": lensInfo[i].getAngularPosition().getY()
            })

        return params

cdef class MassSheetLens(GravitationalLens):
    """Gravitational lens effect of "n infinitely large mass sheet of constant surface mass density"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.MassSheetLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        try:
            # Here we also check for Dd, even though it's not used. This is because
            # it is added in __init__ and we don't want this to be the cause of failure
            GravitationalLens._checkParams(params, ["density", "Dd"])
            return new gravitationallens.MassSheetLensParams(params["density"])
        except LensException as e:
            try:
                GravitationalLens._checkParams(params, ["Dd", "Ds", "Dds"])
                return new gravitationallens.MassSheetLensParams(params["Dd"], params["Ds"], params["Dds"])
            except LensException as e2:
                raise LensException("Parameters must be either 'density', or 'Ds' and 'Dds'")

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have either the following entry:

           * 'density':  the surface mass density of the lens
        
           or the following entries:

           * 'Ds'  : the angular diameter distance between observer and source plane
           * 'Dds' : the angular diameter distance between lens and source plane

           In the second case, the critical density corresponding to ``Dd``, ``Dds`` and ``Ds``
           is used.
        """
        super(MassSheetLens, self).__init__(_gravLensRndId)
        if params is not None:
            params = copy.deepcopy(params)
            params["Dd"] = Dd
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.MassSheetLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.MassSheetLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a MassSheetLens")

        return { "density": pParams.getDensity() }

cdef class CompositeLens(GravitationalLens):
    """This models a gravitational lens that consists of other known lens models"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.CompositeLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef gravitationallens.CompositeLensParams *pLensParams = new gravitationallens.CompositeLensParams()
        cdef gravitationallens.GravitationalLens *pLens = NULL

        try:
            if not isinstance(params, list):
                raise LensException("Composite lens parameters must be a list of individual parameters")

            if len(params) <= 0:
                raise LensException("CompositeLens parameters must be a list")

            for p in params:
                GravitationalLens._checkParams(p, [ "lens", "factor", "angle", "x", "y"])
                pLens = GravitationalLens._getLens(p["lens"])

                if not pLensParams.addLens(p["factor"], Vector2Dd(p["x"], p["y"]), p["angle"], deref(pLens)):
                    raise LensException(S(pLensParams.getErrorString()))
        except Exception as e:
            del pLensParams
            raise

        return pLensParams

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a list of which each entry is a dictionary with the following entries:

           * 'lens':  the gravitational lens object of this component
           * 'factor': the effect of this component can be amplified or reduced by setting 
             this factor. Set to 1 to use the component as-is.
           * 'angle': the angle in degrees that the component should be rotated (counter 
             clockwise)
           * 'x':     the center of the component is now placed at this X-position
           * 'y':     the center of the component is now placed at this Y-position

        **Warning!** 
        No check is done to make sure that the same ``Dd`` parameter from the constructor
        is also used in the individual components. Make sure that this is the case, or
        deal with undefined behaviour!
        """
        super(CompositeLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.CompositeLens *pCompLens
        cdef const gravitationallens.GravitationalLens *pSubLensConst
        cdef unique_ptr[gravitationallens.GravitationalLens] pSubLens
        cdef int i, num
        cdef Vector2Dd pos
        
        self._check()
        pCompLens = dynamic_cast[gravitationallens.CompositeLensPtr](GravitationalLens._getLens(self))
        if pCompLens == NULL:
            raise LensException("Unexpected: this does not seem to be a CompositeLens")

        params = [ ]
        num = pCompLens.getNumberOfSubLenses()
        for i in range(num):
            pSubLensConst = pCompLens.getSubLens(i)
            if pSubLensConst == NULL:
                raise LensException("Unexpected: a sublens is NULL")

            pSubLens = pSubLensConst.createCopy()
            if pSubLens.get() == NULL:
                raise LensException("Unexpected: can't create a copy of a sublens")

            pos = pCompLens.getSubLensPosition(i)
            params.append({
                "lens": GravitationalLens._finalizeLoadedLens(pSubLens),
                "factor": pCompLens.getSubLensFactor(i),
                "angle": pCompLens.getSubLensAngleOriginal(i),
                "x": pos.getX(),
                "y": pos.getY()
            })

        return params

    def findCLSubroutines(self, deflectionScale, potentialScale, bool derivatives, bool potential):
        """findCLSubroutines(deflectionScale, potentialScale, derivatives, potential)

        Analyzes the current CompositeLens instance, returns a tuple of:

         - the recursion level needed (how many other CompositeLens levels 
           it contains)
         - a dictionary where the keys are subroutine names needed by the full
           program (e.g. a plummer lens function) and the values are the actual
           OpenCL programs for these subroutines
         - an array where each entry is either the subroutine name for a sublens
           or `None` if not needed.
        
        The recursion level and array can be fed into the :func:`getCompositeCLProgram`
        function.
        """

        cdef int maxRecursion, i
        cdef cmap[string, string] subRoutineCodes
        cdef cmap[string, string].iterator it
        cdef vector[string] otherRoutineNames
        cdef const gravitationallens.CompositeLens *pLens = gravitationallens.CompositeLens.cast(self._lens())

        if not pLens:
            raise LensException("Internal error: this lens does not seem to be a CompositeLens")
        
        maxRecursion = pLens.findCLSubroutines(deflectionScale, potentialScale, subRoutineCodes, otherRoutineNames, derivatives, potential)

        subCodes = { }
        it = subRoutineCodes.begin()
        while it != subRoutineCodes.end():
            subCodes[S(deref(it).first)] = S(deref(it).second)
            postincrement(it)

        otherRout = []
        for i in range(otherRoutineNames.size()):
            otherRout.append(S(otherRoutineNames[i]) if otherRoutineNames[i].length() > 0 else None)

        return maxRecursion, subCodes, otherRout

    @staticmethod
    def getCompositeCLProgram(deflectionScale, potentialScale, otherRoutineNames, int maxRecursion, bool derivatives, bool potential):
        """getCompositeCLProgram(deflectionScale, potentialScale, otherRoutineNames, maxRecursion, derivatives, potential)

        Returns the name of the OpenCL function as well a the OpenCL program itself for a
        CompositeLens where `maxRecursion` levels are needed, and calls to the subroutine 
        names in `otherRoutineNames` need to be present. These two parameters are returned 
        by the :func:`findCLSubroutines` function. The code for these other subroutines is 
        not included."""

        cdef vector[string] otherNames
        cdef string subName, prog, empty

        for s in otherRoutineNames:
            if s:
                otherNames.push_back(B(s))
            else:
                otherNames.push_back(empty)

        prog = gravitationallens.CompositeLens.getCLProgram_static(deflectionScale, potentialScale, subName, otherNames, maxRecursion, derivatives, potential)
        return S(subName),S(prog)
        

cdef class MassDiskLens(GravitationalLens):
    """This models the gravitational lens effect of a disk with a constant
    density that's centered on the origin."""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.MassDiskLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        try:
            # Here we also check for Dd, even though it's not used. This is because
            # it is added in __init__ and we don't want this to be the cause of failure
            GravitationalLens._checkParams(params, ["density", "radius", "Dd" ]) 
            return new gravitationallens.MassDiskLensParams(params["density"], params["radius"])
        except LensException as e:
            try:
                GravitationalLens._checkParams(params, ["Dd", "Ds", "Dds", "radius"])
                return new gravitationallens.MassDiskLensParams(params["Dd"], params["Ds"], params["Dds"], params["radius"])
            except LensException as e2:
                raise LensException("Parameters must be either 'density' and 'radius', or 'Ds', 'Dds' and 'radius'")

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have either the following entry:

           * 'density': the surface mass density of the lens
           * 'radius':  the angular radius of the disk
        
           or the following entries:

           * 'Ds'  : the angular diameter distance between observer and source plane
           * 'Dds' : the angular diameter distance between lens and source plane
           * 'radius':  the angular radius of the disk

           In the second case, the critical density corresponding to ``Dd``, ``Dds`` and ``Ds``
           is used.
        """
        super(MassDiskLens, self).__init__(_gravLensRndId)
        if params is not None:
            params = copy.deepcopy(params)
            params["Dd"] = Dd
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.MassDiskLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.MassDiskLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a MassDiskLens")

        return { "density": pParams.getDensity(), "radius": pParams.getAngularRadius() }

cdef class PolynomialMassProfileLens(GravitationalLens):
    """With this gravitational lens model, a circularly symmetric lens centered on
    the origin can be modeled based on a specified enclosed mass profile. This
    profile can be specified by different pieces, each one being a polynomial.
    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PolynomialMassProfileLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef gravitationallens.PolynomialMassProfileLensParams *lensParams = new gravitationallens.PolynomialMassProfileLensParams()
        cdef gravitationallens.PolynomialPart *polyPart = NULL

        try:
            for p in params:
                GravitationalLens._checkParams(p, [ "xoffset", "yoffset", "xscale", "yscale", "xend", "coeffs" ])

                del polyPart
                polyPart = new gravitationallens.PolynomialPart(p["xoffset"], p["yoffset"], p["xscale"], p["yscale"], p["xend"], p["coeffs"])
                lensParams.addPolynomialPart(deref(polyPart))
        except:
            del lensParams
            del polyPart
            raise

        del polyPart
        return lensParams

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a list of which each entry is a dictionary with the following entries:

           * 'xoffset': see below
           * 'yoffset': see below
           * 'xscale': see below
           * 'yscale': see below
           * 'xend': see below
           * 'coeffs': see below
            
           The 'xend' values specify the x values (the :math:`\theta` radii) where each polynomial part
           is valid: the current profile part is used if x is between the
           previous list entry's 'xend' value (or 0 if it doesn't exist) and the current
           'xend' value.

           If it is used, the mass profile's value M(x) will be calculated as:
           
           .. code-block:: python

                M(x) = yscale * ( ... + coeffs[k] * ((x-xoffset)/scale)^k + ... ) + yoffset

        """
        super(PolynomialMassProfileLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.PolynomialMassProfileLensParamsPtrConst pParams
        cdef vector[gravitationallens.PolynomialPart] parts
        cdef vector[double] partCoeffs
        cdef int i, l, pl
        
        self._check()
        pParams = dynamic_cast[gravitationallens.PolynomialMassProfileLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a PolynomialMassProfileLens")

        params = [ ]
        parts = pParams.getPolynomialParts()
        l = parts.size()
        for i in range(l):
            coeffs = [ ]
            partCoeffs = parts[i].getCoefficients()
            for j in range(len(partCoeffs)):
                coeffs.append(partCoeffs[j])

            params.append({
                'xoffset': parts[i].getXOffset(),
                'yoffset': parts[i].getYOffset(),
                'xscale': parts[i].getXScale(),
                'yscale': parts[i].getYScale(),
                'xend': parts[i].getEndPosition(),
                'coeffs': coeffs
            })

        return params

cdef class TimeDelayAdjustLens(GravitationalLens):
    """A circularly symmetric lens with which the timedelay in the central
    part can be changed by a specific amount. An example can be found
    in the notebook `timedelayadjust.ipynb <_static/timedelayadjust.ipynb>`_"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PolynomialMassProfileLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, [ "theta1", "theta2", "z", "dt" ])
        return new gravitationallens.PolynomialTimeDelayLensParams(params["theta1"], params["theta2"], params["dt"], params["z"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'theta1': angular radius withing which the lensing potential has a constant value, thereby
             adjusting the time delay.
           * 'theta2': the angular radius beyond which the lensing potential is zero. In between the
             potential will change smoothly.
           * 'dt': the amount with which the time delay in the central region should be modified,
             determining the value of the lensing potential within 'theta1'.
           * 'z': the redshift of the lens.
        """
        super(TimeDelayAdjustLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        raise LensException("The parameters for a TimeDelayAdjustLens cannot be retrieved, as it is actually a PolynomialMassProfileLens. " +
              "Use getPolynomialLensParameters() to obtain the parameters for that type of lens.")

    def getPolynomialLensParameters(self):
        b = self.toBytes()
        l = GravitationalLens.fromBytes(b)
        return l.getLensParameters()

cdef class ZeroMassLens(GravitationalLens):
    """This produces a symmetric gravitational lens with zero total mass,
    as described in the article below.

    **References**

     * `Liesenborgs, J., De Rijcke, S., Dejonghe, H., Bekaert, P.,     
       Non-parametric strong lens inversion of Cl 0024+1654: illustrating the monopole degeneracy,
       Monthly Notices of the Royal Astronomical Society, Volume 389, Issue 1, pp. 415-422, September 2008.
       <http://adsabs.harvard.edu/abs/2008MNRAS.389..415L>`_

    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PolynomialMassProfileLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        GravitationalLens._checkParams(params, [ "Dd", "density", "radius", "zeropoint" ])
        return new gravitationallens.PolynomialZeroMassLensParams(params["Dd"], params["density"], params["radius"], params["zeropoint"])

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have the following entries:

           * 'density': the surface mass density at the center of the lens
           * 'radius': the angular radius beyond which the total mass of the lens will be zero
           * 'zeropoint': a value between 0 and 1 that specified the position (as a fraction of 
             the radius) where the surface mass density becomes zero and the negative density 
             begins.
        """
        super(ZeroMassLens, self).__init__(_gravLensRndId)
        if params is not None:
            params = copy.deepcopy(params)
            params["Dd"] = Dd
        self._lensInit(Dd, params)

    def getLensParameters(self):
        raise LensException("The parameters for a ZeroMassLens cannot be retrieved, as it is actually a PolynomialMassProfileLens. " +
              "Use getPolynomialLensParameters() to obtain the parameters for that type of lens.")

    def getPolynomialLensParameters(self):
        b = self.toBytes()
        l = GravitationalLens.fromBytes(b)
        return l.getLensParameters()

cdef class MassDiskLensSmoothed(GravitationalLens):
    r"""This is the same a a mass disk (:class:`MassDiskLens`), but with a
    smoother edge.
    
    The density profile is:

    .. math::

        \Sigma(\theta) = \left\{\begin{array}{lr}
            \Sigma_0 & \theta < \theta_1 \\
            \Sigma_0 \left[ 2\left(\frac{\theta-\theta_1}{\theta_2-\theta_1}\right)^3
                            -3\left(\frac{\theta-\theta_1}{\theta_2-\theta_1}\right)^2
                            +1
                     \right] & \theta \in [\theta_1, \theta_2] \\
            0 & \theta > \theta_2 
            \end{array}\right.
"""
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PolynomialMassProfileLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef gravitationallens.PolynomialMassProfileLensParams *lensParams = new gravitationallens.PolynomialMassProfileLensParams()
        cdef double dens = 0
        cdef double radius = 0
        cdef double endradius = 0
        cdef double Dd = 0
        cdef double factor = 0
        cdef double x1 = 0
        try:
            # Note that Dd is added in the constructor
            try:
                GravitationalLens._checkParams(params, ["Dd", "density", "radius", "endradius"])
                dens = params["density"]

            except LensException as e:
                try:
                    GravitationalLens._checkParams(params, ["Dd", "Ds", "Dds", "radius", "endradius"])
                except LensException as e:
                    raise LensException("Parameters must be either 'density', 'radius' and 'endradius', or 'Ds', 'Dds', 'radius' and 'endradius'")

                dens = constants.SPEED_C**2/(4.0*np.pi*constants.CONST_G*params["Dd"])*(params["Ds"]/params["Dds"])

            Dd = params["Dd"]
            radius = params["radius"]
            endradius = params["endradius"]

            if radius <= 0 or endradius <= radius:
                raise LensException("The disk radius must be positive and the end radius must be larger")
    
            x1 = radius/(endradius-radius)
            factor = 2.0*np.pi*Dd*Dd

            lensParams.addPolynomialPart(gravitationallens.PolynomialPart(0, 0, 1, factor*dens, radius, [0, 0, 0.5]))
            lensParams.addPolynomialPart(gravitationallens.PolynomialPart(radius, factor*radius**2/2.0, endradius-radius, factor*dens*(endradius-radius)**2, endradius, [ 0, x1, 0.5, -x1, (-0.75+x1/2.0), 0.4 ]))

        except:
            del lensParams
            raise

        return lensParams

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params is a dictionary that should have either the following entry:

           * 'density': the surface mass density of the lens
           * 'radius': the angular radius of the disk itself
           * 'endradius': the angular radius where the density becomes zero
        
           or the following entries:

           * 'Ds'  : the angular diameter distance between observer and source plane
           * 'Dds' : the angular diameter distance between lens and source plane
           * 'radius': the angular radius of the disk itself
           * 'endradius': the angular radius where the density becomes zero

           In the second case, the critical density corresponding to ``Dd``, ``Dds`` and ``Ds``
           is used.
        """
        super(MassDiskLensSmoothed, self).__init__(_gravLensRndId)
        if params is not None: 
            params = copy.deepcopy(params)
            params["Dd"] = Dd
        self._lensInit(Dd, params)

    def getLensParameters(self):
        raise LensException("The parameters for a MassDiskLensSmoothed cannot be retrieved, as it is actually a PolynomialMassProfileLens. " +
              "Use getPolynomialLensParameters() to obtain the parameters for that type of lens.")

    def getPolynomialLensParameters(self):
        b = self.toBytes()
        l = GravitationalLens.fromBytes(b)
        return l.getLensParameters()

cdef class MultipleWendlandLens(GravitationalLens):
    """This creates a lens that's composed of the Wendland basis functions described
    in the LensPerfect article. You can either specify the basis functions yourself,
    or you can specify a list of points at which certain deflection angles should be
    present in the resulting lens. 

    An example notebook can be found here: `wendland.ipynb <_static/wendland.ipynb>`_, the lens file
    used in the notebook is the following: `reallens_nosheet.lensdata <_static/reallens_nosheet.lensdata>`_

    **References**

     * `Coe, D., Fuselier, E., Benítez, N., Broadhurst, T., Frye, B., Ford, H.,
       LensPerfect: Gravitational Lens Mass Map Reconstructions Yielding Exact Reproduction of All Multiple Images,
       The Astrophysical Journal, Volume 681, Issue 2, pp. 814-830, July 2008.
       <http://adsabs.harvard.edu/abs/2008ApJ...681..814C>`_

    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.MultipleWendlandLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:

        cdef gravitationallens.MultipleWendlandLensParams *pParams = NULL
        cdef vector[gravitationallens.WendlandLensInfo] phiXArray;
        cdef vector[gravitationallens.WendlandLensInfo] phiYArray;
        cdef vector[Vector2Dd] points;
        cdef vector[Vector2Dd] angles;
        cdef gravitationallens.WendlandLensInfo w;

        modePhiXPhiY = False
        raiseExcept = False
        try:
            # Here we also check for Dd, even though it's not used. This is because
            # it is added in __init__ and we don't want this to be the cause of failure
            GravitationalLens._checkParams(params, ["phix", "phiy"])
            modePhiXPhiY = True
        except LensException as e:
            pass

        if not modePhiXPhiY:
            try:
                GravitationalLens._checkParams(params, ["points", "angles", "scale"])
            except LensException as e:
                raiseExcept = True

        if raiseExcept:
            raise LensException("Parameters must be either contain arrays 'phix' and 'phiy', or arrays 'points', 'angles' and a value for 'scale'")

        if modePhiXPhiY:
            for p in params["phix"]:
                GravitationalLens._checkParams(p, ["weight", "scale", "position"])
                w = gravitationallens.WendlandLensInfo(p["weight"], p["scale"], Vector2Dd(p["position"][0], p["position"][1]))
                phiXArray.push_back(w)
            for p in params["phiy"]:
                GravitationalLens._checkParams(p, ["weight", "scale", "position"])
                w = gravitationallens.WendlandLensInfo(p["weight"], p["scale"], Vector2Dd(p["position"][0], p["position"][1]))
                phiYArray.push_back(w)

            return new gravitationallens.MultipleWendlandLensParams(phiXArray, phiYArray) 
        else:
            for p in params["points"]:
                points.push_back(Vector2Dd(p[0],p[1]))
            for p in params["angles"]:
                angles.push_back(Vector2Dd(p[0],p[1]))

            pParams = new gravitationallens.MultipleWendlandLensParams()
            if not pParams.matchDeflections(points, angles, params["scale"]):
                errStr = S(pParams.getErrorString())
                del pParams
                raise LensException(errStr)

            return pParams

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: should be one of two cases:

            1. a dictionary with entries 'phix' and 'phiy'. Each of those should be a list
               of dictionaries, each describing a Wendland basis function that's oriented
               along the x- or y-axis. These dictionaries should have the entries:

                - 'weight': a weight for the basis function
                - 'scale': an angular scale for the basis function
                - 'position': the angular position of the basis function

            2. a dictionary containing entries

                - 'points': a list of angular positions at which the deflection angle will
                  be specified
                - 'angles': a list of equal length as the 'points' list, describing the
                  deflection angle the lens should have at that point.
                - 'scale': an angular scale of the basis functions

                If this method is used, the basis functions will be determined automatically
                in such a way that they generate the specified deflection angles at the specified
                points, where each basis function has the angular scale given by 'scale'.
        """
        super(MultipleWendlandLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.MultipleWendlandLensParamsPtrConst pParams
        cdef vector[gravitationallens.WendlandLensInfo] lensInfo
        cdef int i, l
        
        self._check()
        pParams = dynamic_cast[gravitationallens.MultipleWendlandLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a MultipleWendlandLens")

        phiX = [ ]
        phiY = [ ]

        lensInfo = pParams.getPhiXInfo()
        l = lensInfo.size()
        for i in range(l):
            phiX.append({
                "weight": lensInfo[i].getHeightFactor(),
                "scale": lensInfo[i].getAngularScale(),
                "position": [ lensInfo[i].getAngularPosition().getX(),
                              lensInfo[i].getAngularPosition().getY() ]
            })

        lensInfo = pParams.getPhiYInfo()
        l = lensInfo.size()
        for i in range(l):
            phiY.append({
                "weight": lensInfo[i].getHeightFactor(),
                "scale": lensInfo[i].getAngularScale(),
                "position": [ lensInfo[i].getAngularPosition().getX(),
                              lensInfo[i].getAngularPosition().getY() ]
            })

        return { "phix": phiX, "phiy": phiY }

cdef class DeflectionGridLens(GravitationalLens):
    """Creates a lens based on provided deflection angles, provided on a regular grid.
    When a deflection angle is later retrieved using e.g. 
    :py:meth:`getAlphaVector<grale.lenses.GravitationalLens.getAlphaVector>`,
    each component will be interpolated bilinearly based on the provided deflections. This
    type of lens is what is created when 
    :py:meth:`LensPlane.createDeflectionGridLens<grale.images.LensPlane.createDeflectionGridLens>`
    is called.

    Note that this is only an approximation and should be used with some limitations
    in mind. For one thing, it will not have a valid lens potential associated with it.
    Also, as the deflection angles rely on interpolation, they will in general *not* be
    curl-free.

    An example notebook can be found here: `deflectiongridlens.ipynb <_static/deflectiongridlens.ipynb>`_, 
    the lens file used in the notebook is the following: `reallens_nosheet.lensdata <_static/reallens_nosheet.lensdata>`_
    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.DeflectionGridLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:

        cdef Vector2Dd bl, tr;
        cdef np.ndarray[double,ndim=3] deflectionAngles
        cdef vector[double] alphaX, alphaY
        cdef int x, y, w, h;

        GravitationalLens._checkParams(params, ["bottomleft", "topright", "angles"])
        tr = Vector2Dd(params["topright"][0], params["topright"][1])
        bl = Vector2Dd(params["bottomleft"][0], params["bottomleft"][1])
        deflectionAngles = params["angles"]

        w = deflectionAngles.shape[1]
        h = deflectionAngles.shape[0]
        for y in range(h):
            for x in range(w):
                alphaX.push_back(deflectionAngles[y,x,0])
                alphaY.push_back(deflectionAngles[y,x,1])

        return new gravitationallens.DeflectionGridLensParams(alphaX, alphaY, w, h, bl, tr)

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

             - 'angles': a NumY*NumX*2 numpy array containing the deflection angles.
             - 'bottomleft': the bottom-left coordinate, which will have the deflection
               angle stored at position (0, 0) in the 'angles' grid.
             - 'topright': the top-right coordinat, which will have the deflection angle
               stored at position (NumY-1, NumX-1) in the 'angles' grid.
        """
        super(DeflectionGridLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.DeflectionGridLensParamsPtrConst pParams
        cdef vector[double] alphaX
        cdef vector[double] alphaY
        cdef int w, h, x, y
        cdef np.ndarray[double,ndim=3] deflectionAngles
        
        self._check()
        pParams = dynamic_cast[gravitationallens.DeflectionGridLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a DeflectionGridLens")

        w = pParams.getWidth()
        h = pParams.getHeight()
        deflectionAngles = np.zeros([h, w, 2], dtype = np.double)
        alphaX = pParams.getAlphaX()
        alphaY = pParams.getAlphaY()

        if <int>alphaX.size() != w*h or alphaX.size() != alphaY.size():
            raise LensException("Internal error: alphaX and alphaY sizes are not compatible with grid size")

        for y in range(h):
            for x in range(w):
                deflectionAngles[y,x,0] = alphaX[x+y*w]
                deflectionAngles[y,x,1] = alphaY[x+y*w]

        return { 
            "angles": deflectionAngles, 
            "bottomleft": [ pParams.getBottomLeft().getX(), pParams.getBottomLeft().getY() ],
            "topright": [ pParams.getTopRight().getX(), pParams.getTopRight().getY() ]
        }

cdef class NFWLens(GravitationalLens):
    r"""A lens model for a 2D symmetric lens profile, based on a 3D symmetric 
Navarro-Frenk-White (NFW) mass distribution. 

We'll be using the following definitions of the helper functions :math:`F(x)`,
:math:`G(x)`, :math:`H(x)`, :math:`A(x)` and :math:`B(x)` as follows

.. math::

    F(x) = \left\{
    \begin{array}{ll} 
     	    \frac{1}{\sqrt{1-x^2}}\textrm{atanh}{\sqrt{1-x^2}} & x < 1 \\
   	        1 & x = 1 \\
   	        \frac{1}{\sqrt{x^2-1}}\textrm{atan}{\sqrt{x^2-1}} & x > 1 \\
   	\end{array} \right.
    \Rightarrow \frac{d F}{d x} = \frac{1-x^2 F(x)}{x(x^2-1)}

.. math::

    G(x) = \frac{1-F(x)}{x^2-1} \qquad 
    H(x) = \textrm{ln}\left(\frac{x}{2}\right) + F(x) \Rightarrow \frac{d H}{d x} = G(x) x

.. math::
    
    B(x) = \left\{
       \begin{array}{ll} 
   	        \textrm{ln}^2\left(\frac{x}{2}\right) - \textrm{atanh}^2\sqrt{1-x^2} & x < 1 \\
   	        \textrm{ln}^2 2 & x = 1 \\
   	        \textrm{ln}^2\left(\frac{x}{2}\right) + \textrm{atan}^2\sqrt{x^2-1} & x > 1 \\
   	    \end{array} \right.
    \qquad 
    A(x) = \frac{d H}{d x} 
    \Rightarrow \frac{1}{2} \frac{d B}{d x} = A(x)

The following relations then hold for the 2D surface mass density :math:`\Sigma(\theta)` and the 2D
integrated mass profile :math:`\textrm{M}(\theta)`

.. math::

    \Sigma(\theta) = 2 r_s \rho_s G\left(\frac{\theta}{\theta_s}\right)
    \qquad
    \textrm{M}(\theta) = 4\pi r_s^3\rho_s H\left(\frac{\theta}{\theta_s}\right)

Being a symmetric lens, this :math:`\textrm{M}(\theta)` then leads to the deflection angles

.. math::

    \vec{\hat{\alpha}}\left(\vec{\theta}\right) = \frac{16 \pi G r_s^2 \rho_s}{c^2} A\left(\frac{\theta}{\theta_s}\right) \frac{\vec{\theta}}{\theta}

while the (unscaled) lensing potential becomes

.. math::

    \psi(\theta) = \frac{D_{ds}}{D_s} \frac{8 \pi G r_s^2 \rho_s \theta_s}{c^2} B\left(\frac{\theta}{\theta_s}\right)

**References**

 * `Keeton, C., A Catalog of Mass Models for Gravitational Lensing, February 2001.
   <http://adsabs.harvard.edu/abs/2001astro.ph..2341K>`_
 * `Wright, C., Brainerd, T., Gravitational Lensing by NFW Halos, The Astrophysical Journal, Volume 534, Issue 1, pp. 34-40, May 2000.
   <http://adsabs.harvard.edu/abs/2000ApJ...534...34W>`_


    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.NFWLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double rho_s, theta_s

        GravitationalLens._checkParams(params, ["rho_s", "theta_s"])
        rho_s = params["rho_s"]
        theta_s = params["theta_s"]
        return new gravitationallens.NFWLensParams(rho_s, theta_s)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

            - 'rho_s': a density scale for the NFW profile, corresponding to :math:`\rho_s` in the
              equations above
            - 'theta_s': an angular scale factor for the NFW profile, corresponding to
              :math:`\theta_s` in the equations above
        """
        super(NFWLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.NFWLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.NFWLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a NFWLens")

        return { "rho_s": pParams.get3DDensityScale(), "theta_s": pParams.getAngularRadiusScale() }

cdef class EllipticNFWLens(GravitationalLens):
    r"""An elliptic generalization of the :class:`NFWLens` lens model. The relations
    between the circularly symmetric profile and the elliptic generalization are
    described in the reference.

    The degree of ellipticity is encoded by a parameter :math:`q`, so that for the
    2D mass density the relation holds:

    .. math::

        \Sigma(\theta_x, \theta_y) = \Sigma_{\rm circular}\left(\theta_x, \frac{\theta_y}{q}\right)

    **References**

     * `Keeton, C., A Catalog of Mass Models for Gravitational Lensing, February 2001.
       <http://adsabs.harvard.edu/abs/2001astro.ph..2341K>`_
    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.EllipticNFWLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double rho_s, theta_s, q

        GravitationalLens._checkParams(params, ["rho_s", "theta_s", "q"])
        rho_s = params["rho_s"]
        theta_s = params["theta_s"]
        q = params["q"]
        return new gravitationallens.EllipticNFWLensParams(rho_s, theta_s, q)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

            - 'rho_s': a density scale for the NFW profile, corresponding to :math:`\rho_s` in the
              equations in :class:`NFWLens`
            - 'theta_s': an angular scale factor for the NFW profile, corresponding to
              :math:`\theta_s` in the equations in :class:`NFWLens`
            - 'q': a measure of the ellipticity of the lens, where 1 describes again the
              circular profile.
        """
        super(EllipticNFWLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.EllipticNFWLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.EllipticNFWLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of an EllipticNFWLens")

        return { 
            "rho_s": pParams.get3DDensityScale(), 
            "theta_s": pParams.getAngularRadiusScale(),
            "q": pParams.getEllipticity()
            }

cdef class SersicLens(GravitationalLens):
    r"""A lens model with a 2D projected mass density based on the Sérsic profile. This
    mass density is the following:

    .. math::

        \Sigma(\theta) = \Sigma_{\rm central}\exp\left[\left(\frac{\theta}{\theta_s}\right)^{\frac{1}{n}}\right]
    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.SersicLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double sigma0, thetaScale, index

        GravitationalLens._checkParams(params, ["centraldensity", "scale", "index"])
        sigma0 = params["centraldensity"]
        thetaScale = params["scale"]
        index = params["index"]

        return new gravitationallens.SersicLensParams(sigma0, thetaScale, index)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

            - 'centraldensity': the value of the 2D density at the center of the
              lens corresponding to :math:`\Sigma_{\rm central}` in the equation above
            - 'scale': an angular scale of the lens, corresponding to :math:`\theta_s` 
              in the equation above
            - 'index': the so-called Sérsic index, corresponding to :math:`n` in the
              equation above
        """
        super(SersicLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.SersicLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.SersicLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a SersicLens")

        return { 
            "centraldensity": pParams.getCentralDensity(), 
            "scale": pParams.getAngularScale(),
            "index": pParams.getSersicIndex()
        }

cdef class EllipticSersicLens(GravitationalLens):
    r"""An elliptic generalization of the :class:`SersicLens` lens model. The relations
    between the circularly symmetric profile and the elliptic generalization are
    described in the reference.

    The degree of ellipticity is encoded by a parameter :math:`q`, so that for the
    2D mass density the relation holds:

    .. math::

        \Sigma(\theta_x, \theta_y) = \Sigma_{\rm circular}\left(\theta_x, \frac{\theta_y}{q}\right)

    **References**

     * `Keeton, C., A Catalog of Mass Models for Gravitational Lensing, February 2001.
       <http://adsabs.harvard.edu/abs/2001astro.ph..2341K>`_
    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.EllipticSersicLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double sigma0, thetaScale, index, q

        GravitationalLens._checkParams(params, ["centraldensity", "scale", "index", "q"])
        sigma0 = params["centraldensity"]
        thetaScale = params["scale"]
        index = params["index"]
        q = params["q"]

        return new gravitationallens.EllipticSersicLensParams(sigma0, thetaScale, index, q)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

            - 'centraldensity': the value of the 2D density at the center of the
              lens corresponding to :math:`\Sigma_{\rm central}` in the equation 
              in :class:`SersicLens`
            - 'scale': an angular scale of the lens, corresponding to :math:`\theta_s` 
              in the equation in :class:`SersicLens`
            - 'index': the so-called Sérsic index, corresponding to :math:`n` in the
              equation in :class:`SersicLens`
            - 'q': a measure of the ellipticity of the lens, where 1 describes again the
              circular profile.
        """
        super(EllipticSersicLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.EllipticSersicLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.EllipticSersicLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of an EllipticSersicLens")

        return { 
            "centraldensity": pParams.getCentralDensity(), 
            "scale": pParams.getAngularScale(),
            "index": pParams.getSersicIndex(),
            "q": pParams.getEllipticity()
        }

cdef class ProfileLens(GravitationalLens):
    """A circularly symmetric lens based on a discrete density profile. In
    between the points of the profile, the density is linearly interpolated.

    An example notebook can be found here: `profilelens.ipynb <_static/profilelens.ipynb>`_.
    """
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.ProfileLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double radius
        cdef vector[double] profile;

        GravitationalLens._checkParams(params, ["radius", "profile"])
        radius = params["radius"]
        for x in params["profile"]:
            profile.push_back(x)

        return new gravitationallens.ProfileLensParams(radius, profile)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

            - 'radius': the profile will contain density values for evenly spaced points between
              zero and this end radius. Beyond this radius, the model does not contain any mass.
            - 'profile': a list of density values, the first of which corresponds to the central
              density of the circularly symmetric lens, and the last corresponds to the density
              at the specified end radius.
        """
        super(ProfileLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.ProfileLensParamsPtrConst pParams
        cdef vector[double] profileVector
        cdef int i, l
        
        self._check()
        pParams = dynamic_cast[gravitationallens.ProfileLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a ProfileLens")

        profile = [ ]
        profileVector = pParams.getProfile()
        l = profileVector.size()
        for i in range(l):
            profile.append(profileVector[i])

        return { "radius": pParams.getEndRadius(), "profile": profile }

cdef class PIEMDLens(GravitationalLens):
    r"""A lens with a PIEMD profile with core and cut radius. Although this is
    usually referred to as a PIEMD, it is not a PIEMD as described in the article
    of Kassiola and Kovner. In Eliasdottir et al, is is called a dPIE to make the
    distinction more clear.

    **References**

     * `Eliasdottir, A., et al, Where is the matter in the Merging Cluster Abell 2218? <http://adsabs.harvard.edu/abs/2007arXiv0710.5636E>`_
     * `Kassiola, A., Kovner, I., Analytic lenses with elliptic mass distributions, pseudo-isothermal and others, instead of elliptic potentials <http://adsabs.harvard.edu/abs/1993LIACo..31..571K>`_
    """

    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PIEMDLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double sigma0, coreRadius, scaleRadius, epsilon

        GravitationalLens._checkParams(params, ["centraldensity", "coreradius", "scaleradius", "epsilon"])
        sigma0 = params["centraldensity"]
        coreRadius = params["coreradius"]
        scaleRadius = params["scaleradius"]
        epsilon = params["epsilon"]

        return new gravitationallens.PIEMDLensParams(sigma0, coreRadius, scaleRadius, epsilon)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

           * 'centraldensity': the central projected density :math:`\Sigma_0`
           * 'coreradius': the core radius :math:`a` of the profile
           * 'scaleradius': the scale radius :math:`s` of the profile
           * 'epsilon': a measure for the ellipticity of the lens

        """
        super(PIEMDLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.PIEMDLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.PIEMDLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a PIEMDLens")

        return { 
            "centraldensity": pParams.getCentralDensity(), 
            "coreradius": pParams.getCoreRadius(),
            "scaleradius": pParams.getScaleRadius(),
            "epsilon": pParams.getEpsilon()
        }

cdef class PIMDLens(GravitationalLens):
    r"""The circularly symmetric version of the PIEMD lens (:class:`PIEMDLens`),
    corresponding to en 'epsilon' parameter of zero.
    """

    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PIMDLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double sigma0, coreRadius, scaleRadius

        GravitationalLens._checkParams(params, ["centraldensity", "coreradius", "scaleradius"])
        sigma0 = params["centraldensity"]
        coreRadius = params["coreradius"]
        scaleRadius = params["scaleradius"]

        return new gravitationallens.PIMDLensParams(sigma0, coreRadius, scaleRadius)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries: 

           * 'centraldensity': the central projected density :math:`\Sigma_0`
           * 'coreradius': the core radius :math:`a` of the profile
           * 'scaleradius': the scale radius :math:`s` of the profile

        """
        super(PIMDLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.PIMDLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.PIMDLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a PIMDLens")

        return { 
            "centraldensity": pParams.getCentralDensity(), 
            "coreradius": pParams.getCoreRadius(),
            "scaleradius": pParams.getScaleRadius()
        }

cdef class AlphaPotLens(GravitationalLens):
    r"""This lens is based on the 'alphapot' softened power law lens model 
    from lenstool/gravlens. The unscaled (:math:`D_{ds}/D_s = 1`) lensing 
    potential is taken to be

    .. math::

        \psi(\vec{\theta}) = b\left(s^2 + \theta_x^2 + \frac{\theta_y^2}{q^2} + K^2\theta_x\theta_y\right)^{\frac{\alpha}{2}}

    Note that because :math:`\kappa(\vec{\theta}) = \frac{1}{2}\left(\frac{\partial^2 \psi}{\partial \theta_x^2} + \frac{\partial^2 \psi}{\partial \theta_y^2}\right)`
    you'll need to do an additional scaling (e.g. by adjusting :math:`b` or using a :class:`CompositeLens`)
    by :math:`(1 \textrm{ arcsec})^{\alpha-2}` to get the same mass density as lenstool (there,
    1 arcsec is used as angular unit).

    **References**

     *  `Keeton, C., gravlens 1.06: Software for Gravitational Lensing, January 2004. <http://www.physics.rutgers.edu/~keeton/gravlens/manual.pdf>`_
    """

    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.AlphaPotLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double b, s, q, K2, alpha

        GravitationalLens._checkParams(params, ["b", "s", "q", "K2", "alpha"])
        b = params["b"]
        s = params["s"]
        q = params["q"]
        K2 = params["K2"]
        alpha = params["alpha"]

        return new gravitationallens.AlphaPotLensParams(b, s, q, K2, alpha)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

           * 'b': the value of :math:`b` in lensing potential
           * 's': the value of :math:`s` in lensing potential
           * 'q': the value of :math:`q` in lensing potential
           * 'K2': the value of :math:`K^2` in lensing potential
           * 'alpha': the value of :math:`\alpha` in the lensing potential

        """
        super(AlphaPotLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.AlphaPotLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.AlphaPotLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a AlphaPotLens")

        return { 
            "b": pParams.getB(), 
            "s": pParams.getS(),
            "q": pParams.getQ(),
            "K2": pParams.getK2(),
            "alpha": pParams.getAlpha()
        }

cdef class HarmonicLens(GravitationalLens):
    r"""Harmonic 2D mass density, similar to what's used in the JPEG article.

    .. math::

        \Sigma(\vec{\theta}) = \Sigma_0 \cos\left(\frac{k}{2}\theta_x + \phi_x\right) \cos\left(\frac{l}{2}\theta_y + \phi_y\right)

    **References:**
     *  `Lam, D., A New Approach to Free-Form Cluster Lens Modeling Inspired by the JPEG Image Compression Method. <https://arxiv.org/abs/1906.00006>`_
    """

    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.HarmonicLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef double sigma0, phiX, phiY, k, l

        GravitationalLens._checkParams(params, ["sigma0", "k", "l", "phi_x", "phi_y"])
        sigma0 = params["sigma0"]
        k = params["k"]
        l = params["l"]
        phiX = params["phi_x"]
        phiY = params["phi_y"]

        return new gravitationallens.HarmonicLensParams(sigma0, k, l, phiX, phiY)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary containing the following entries:

           * 'sigma0': Maximum density :math:`\Sigma_0`
           * 'k': value of 'k' in the formula
           * 'l': value of 'l' in the formula
           * 'phi_x': value of :math:`\phi_x` in the formula
           * 'phi_y': value of :math:`\phi_y` in the formula
        """
        super(HarmonicLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.HarmonicLensParamsPtrConst pParams
        
        self._check()
        pParams = dynamic_cast[gravitationallens.HarmonicLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a HarmonicLens")

        return { 
            "sigma0": pParams.getDensityScale(),
            "k": pParams.getK(),
            "l": pParams.getL(),
            "phi_x": pParams.getPhiX(),
            "phi_y": pParams.getPhiY()
        }

cdef class PotentialGridLensBase:
    cdef unique_ptr[gravitationallens.PotentialGridLensBase] m_pLens
    cdef int m_numX, m_numY
    cdef int m_init

    def __init__(self, double Dd, bottomLeft, topRight, int numX, int numY, values = None):

        if numX < 2 or numY < 2:
            raise LensException("Invalid dimensions specified")
        
        self.m_pLens = unique_ptr[gravitationallens.PotentialGridLensBase](new gravitationallens.PotentialGridLensBase(Dd, Vector2Dd(bottomLeft[0], bottomLeft[1]), Vector2Dd(topRight[0], topRight[1]), numX, numY))
        self.m_numX = numX
        self.m_numY = numY
        self.m_init = 0

        if values is not None:
            self.setValues(values)

    cdef gravitationallens.PotentialGridLensBase * _lens(self):
        return self.m_pLens.get()

    cdef _check(self):
        if not self.m_init:
            raise LensException("No initial potential values were set")

    def getValues(self):
        cdef np.ndarray[double,ndim=1] values
        cdef int num = self._lens().values().size()

        self._check()
                
        if num != self.m_numX*self.m_numY:
            raise LensException("Unexpected: incompatible dimensions")

        values = np.empty([num], dtype = np.double)

        for i in range(num):
            values[i] = self._lens().values()[i]
        
        return values.reshape(self.m_numY, self.m_numX)

    def setValues(self, np.ndarray[double, ndim=2] v):
        cdef np.ndarray[double,ndim=1] linArray
        cdef int num = self.m_numX * self.m_numY
        cdef errut.bool_t r

        if v is None:
            raise LensException("No values specified")

        if not self.m_init:
            self._lens().values().resize(num, 0)
            self.m_init = 1

        if num != self._lens().values().size():
            raise LensException("Unexpected: incompatible dimensions")

        if v.shape[0] != self.m_numY or v.shape[1] != self.m_numX:
            raise LensException("Specified values array has wrong shape")

        linArray = v.reshape((num,))
        for i in range(num):
            self._lens().values()[i] = linArray[i]

        r = self._lens().init()
        if not r.success():
            raise LensException(S(r.getErrorString()))

    # TODO: this is some copy/paste work from GravitationalLens, merge this in a general function
    #       the 'self._check()' could be done as a lambda
    cdef _reshapeAndCall1D(self, functionName, thetas, int coreNumIn, int coreNumOut):
        cdef int totalElements, l, i

        self._check()

        l = len(thetas.shape)
        if l == 1:
            return functionName(thetas)

        if l > 1 and thetas.shape[l-1] == coreNumIn:
            totalElements = 1
            for i in range(l):
                totalElements *= thetas.shape[i]

            outShape = thetas.shape[:-1] 
            if coreNumOut > 1:
                outShape += (coreNumOut,)
            return np.reshape(functionName(np.reshape(thetas,[totalElements])), outShape)

        raise LensException("Bad array dimensions")

    # TODO: same
    cdef _getAlphaVector1(self, np.ndarray[double, ndim=1] thetas):
        
        cdef Vector2Dd alpha
        cdef np.ndarray[double,ndim=1] alphas
        cdef int l,i
        cdef errut.bool_t r

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            alphas = np.zeros([l], dtype = np.double)
            for i in range(0,l,2):

                r = self._lens().getAlphaVector(Vector2Dd(thetas[i], thetas[i+1]), cython.address(alpha))
                if not r.success():
                    raise LensException(S(r.getErrorString()))
                alphas[i] = alpha.getX()
                alphas[i+1] = alpha.getY()
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return alphas

    def getAlphaVector(self, thetas):
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x: self._getAlphaVector1(x), thetas, 2, 2)

    cdef _getAlphaVectorDerivatives1(self, np.ndarray[double, ndim=1] thetas):
        
        cdef double axx = 0, ayy = 0, axy = 0
        cdef np.ndarray[double,ndim=1] derivatives
        cdef int l,i,j
        cdef errut.bool_t r

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            derivatives = np.zeros([(l//2)*3], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                r = self._lens().getAlphaVectorDerivatives(Vector2Dd(thetas[i], thetas[i+1]), axx, ayy, axy)
                if not r.success():
                    raise LensException(S(r.getErrorString()))
                derivatives[j] = axx
                derivatives[j+1] = ayy
                derivatives[j+2] = axy
                j += 3
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return derivatives

    def getAlphaVectorDerivatives(self, thetas):
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x: self._getAlphaVectorDerivatives1(x), thetas, 2, 3)

    cdef _getSurfaceMassDensity1(self, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] dens
        cdef int l,i,j
        cdef errut.bool_t r
        cdef double value = 0

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            dens = np.zeros([(l//2)], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                r = self._lens().getSurfaceMassDensity(Vector2Dd(thetas[i], thetas[i+1]), value)
                if not r.success():
                    raise LensException(S(r.getErrorString()))
                dens[j] = value
                j += 1
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return dens

    def getSurfaceMassDensity(self, thetas):
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x: self._getSurfaceMassDensity1(x), thetas, 2, 1)

    cdef _getProjectedPotential1(self, np.ndarray[double, ndim=1] thetas):
        
        cdef np.ndarray[double,ndim=1] potential
        cdef int l,i,j
        cdef double value = 0
        cdef errut.bool_t r

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            potential = np.zeros([(l//2)], dtype = np.double)
            j = 0
            for i in range(0,l,2):
                r = self._lens().getProjectedPotential(Vector2Dd(thetas[i], thetas[i+1]), cython.address(value))
                if not r.success():
                    raise LensException(S(r.getErrorString()))

                potential[j] = value
                j += 1
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return potential

    def getProjectedPotential(self, thetas):
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x: self._getProjectedPotential1(x), thetas, 2, 1)

cdef class PotentialGridLens(GravitationalLens):
    """Create a lens based on the values of the projected potential (for :math:`D_{ds}/D_s = 1`)
    defined on a grid. In between the grid points, bicubic interpolation is used using
    routines from `GSL <https://www.gnu.org/software/gsl/doc/html/interp.html#id1>`_.
    """

    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.PotentialGridLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef int numX, numY, i, x, y
        cdef np.ndarray[double,ndim=2] values
        cdef vector[double] valueArray
        cdef Vector2Dd topRight, bottomLeft

        GravitationalLens._checkParams(params, ["values", "topright", "bottomleft"])

        values = params["values"]
        numY, numX = values.shape[0], values.shape[1]
        valueArray.resize(numX*numY)

        i = 0
        for y in range(numY):
            for x in range(numX):
                valueArray[i] = values[y,x]
                i += 1
        
        topRight = Vector2Dd(params["topright"][0], params["topright"][1])
        bottomLeft = Vector2Dd(params["bottomleft"][0], params["bottomleft"][1])

        return new gravitationallens.PotentialGridLensParams(bottomLeft, topRight, valueArray, numX, numY)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary consisting of the following entries

           * 'bottomleft': the coordinates of the bottom-left corner of the grid
           * 'topright': the top-right corner of the grid
           * 'values': 2D numpy array containing the potential values, sampled at the grid 
             points. If the shape is `(numY, numX)` (first index describes y-direction)
             then point `[0, 0]` corresponds to the bottom-left corner, and `[numY-1, numX-1]`
             to the top-right corner.
        """
        super(PotentialGridLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.PotentialGridLensParamsPtrConst pParams
        cdef np.ndarray[double,ndim=1] values
        cdef int numX, numY, i, n;
        cdef Vector2Dd bottomLeft, topRight

        self._check()
        pParams = dynamic_cast[gravitationallens.PotentialGridLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a PotentialGridLens")

        numX = pParams.getNumX()
        numY = pParams.getNumY()
        bottomLeft = pParams.getBottomLeft()
        topRight = pParams.getTopRight()

        values = np.empty([numX*numY], dtype = np.double)
        n = pParams.getValues().size()
        if n != numX*numY:
            raise LensException("Unexpected: number of potential values is not equal to numX*numY")

        for i in range(n):
            values[i] = pParams.getValues()[i]

        return {
            "topright": [ topRight.getX(), topRight.getY() ],
            "bottomleft": [ bottomLeft.getX(), bottomLeft.getY() ],
            "values": values.reshape(numY, numX)
        }

cdef class CircularPiecesLens(GravitationalLens):
    r"""Create a lens based on the projected potentials of existing lenses. These 
    potentials will be used exactly in ring-shaped regions (possibly with a
    different offset, and rescaled), and interpolated in between. If :math:`\psi_i`
    is the projected potential for one such component, the scaled/offset version
    is

    .. math::

        \psi_i^s = \mathrm{scale}\times(\psi_i + \mathrm{offset})

    In between such regions, the potentials are merged in the following way:

    .. math::

        \psi_i^{\mathrm{interp}}(\vec{\theta}) = f(x) \psi_i^s(\vec{\theta}) + (1-f(x)) \psi_{i+1}^s(\vec{\theta})

    where

    .. math::

        x = \frac{|\vec{\theta}| - \theta_i^{\rm end}}{\theta_{i+1}^{\rm start} - \theta_i^{\rm end}}

    and :math:`\theta_i^{\rm end}` is the outer radius for the ring-like region of
    potential :math:`i`, and :math:`\theta_{i+1}^{\rm start}` is the inner radius of
    the next ring-like region.
    """

    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.CircularPiecesLens()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef vector[gravitationallens.CircularPieceInfo] pieces
        cdef vector[double] interpolCoeffs
        cdef gravitationallens.GravitationalLens *pLens = NULL
        cdef shared_ptr[gravitationallens.GravitationalLens] lens;

        params = copy.deepcopy(params)
        if isinstance(params, dict) and not "coeffs" in params: # Add a default
            params["coeffs"] = [ 1, 0, 0, -10, 15, -6 ]

        GravitationalLens._checkParams(params, ["pieces", "coeffs"])

        for x in params["coeffs"]:
            interpolCoeffs.push_back(x)

        piecesParams = params["pieces"]
        if not isinstance(piecesParams, list):
            raise LensException("Circular pieces lens parameters must be a list of individual parameters")

        if len(piecesParams) <= 0:
            raise LensException("Circular pieces lens parameters must be a non-empty list")

        for d in piecesParams:
            d = copy.copy(d)
            if not "r0" in d: d["r0"] = 0.0
            if not "r1" in d: d["r1"] = float("inf")
            if not "potentialScale" in d: d["potentialScale"] = 1.0
            if not "potentialOffset" in d: d["potentialOffset"] = 0.0

            allowedKeys = [ "lens", "r0", "r1", "potentialOffset", "potentialScale" ]
            for k in d:
                if not k in allowedKeys:
                    raise LensException("Encountered unknown key '{}', expecting one of {}".format(k, allowedKeys))

            pLens = GravitationalLens._getLens(d["lens"])
            lens.reset(pLens.createCopy().release())

            pieces.push_back(gravitationallens.CircularPieceInfo(lens, d["r0"], d["r1"], d["potentialScale"], d["potentialOffset"]))

        return new gravitationallens.CircularPiecesLensParams(pieces, interpolCoeffs)

    def __init__(self, Dd, params):
        r"""__init__(Dd, params)

        Parameters:
         - Dd is the angular diameter distance to the lens.
         - params: a dictionary with two entries

           * 'coeffs': (optional) a list of numbers that represent the coefficients
             of the interpolation function 
             
                .. math ::

                    f(x) = \sum_k a_k x^k

             The default is ``[ 1, 0, 0, -10, 15, -6  ]``, similar to a 
             `smoothstep-like <https://en.wikipedia.org/wiki/Smoothstep>`_ function.

           * 'pieces': a list of dictionaries with the following entries, representing
             the existing lens potentials and their ring-like regions:
          
             * 'lens': lens model
             * 'r0': start radius (zero if not present)
             * 'r1': end radius (infinity if not present)
             * 'potentialOffset': offset to the lens potential for this lens (zero if not present)
             * 'potentialScale': scale factor for the lens potential for this lens (one if not present)
        """

        super(CircularPiecesLens, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.CircularPiecesLensParamsPtrConst pParams
        cdef const gravitationallens.GravitationalLens *pSubLensConst
        cdef unique_ptr[gravitationallens.GravitationalLens] pSubLens
        
        self._check()
        pParams = dynamic_cast[gravitationallens.CircularPiecesLensParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a CircularPiecesLens")

        piecesParams = [ ]
        for i in range(pParams.getPiecesInfo().size()):
            pSubLensConst = pParams.getPiecesInfo()[i].getLens().get()
            if pSubLensConst == NULL:
                raise LensException("Unexpected: a piece lens is NULL")
            
            pSubLens = pSubLensConst.createCopy()
            if pSubLens.get() == NULL:
                raise LensException("Unextected: can't create a copy of a piece lens")

            piecesParams.append({
                "lens": GravitationalLens._finalizeLoadedLens(pSubLens),
                "r0": pParams.getPiecesInfo()[i].getStartRadius(),
                "r1": pParams.getPiecesInfo()[i].getEndRadius(),
                "potentialOffset": pParams.getPiecesInfo()[i].getPotentialOffset(),
                "potentialScale": pParams.getPiecesInfo()[i].getPotentialScale()
            })
        
        coeffs = []
        for i in range(pParams.getInterpolationFunctionCoefficients().size()):
            coeffs.append(pParams.getInterpolationFunctionCoefficients()[i])

        return { "pieces": piecesParams, "coeffs": coeffs }

cdef class MultiPlaneContainer(GravitationalLens):
    """This 'lens' can be used as a container for multiple lenses as different redshifts"""

    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return new gravitationallens.MultiPlaneContainer()

    cdef gravitationallens.GravitationalLensParams* _allocParams(self, params) except NULL:
        cdef gravitationallens.MultiPlaneContainerParams *pParams = new gravitationallens.MultiPlaneContainerParams()
        cdef gravitationallens.GravitationalLens *pLens = NULL

        try:
            if not isinstance(params, list):
                raise LensException("Multi-plane container parameters must be a list of (lens, redshift) pairs")

            if len(params) <= 0:
                raise LensException("The parameter list does not seem to contain any entries")

            for p in params:
                GravitationalLens._checkParams(p, [ "lens", "z" ])
                pLens = GravitationalLens._getLens(p["lens"])

                if not pParams.add(pLens, p["z"]):
                    raise LensException(S(pParams.getErrorString()))
        except Exception as e:
            del pParams
            raise

        return pParams

    def __init__(self, Dd, params):
        """__init__(Dd, params)

        Parameters:
         - Dd is not used in this container, must be set to 0
         - params is a list of which each entry is a dictionary with the following entries:

           * 'lens':  the gravitational lens object at this redshift
           * 'z':     the redshift for this lens

        **Warning!** 
        No check is done to make sure that the lens distance parameter matches with the
        specified redshift (there's no notion of a cosmological model here)
        """
        super(MultiPlaneContainer, self).__init__(_gravLensRndId)
        self._lensInit(Dd, params)

    def getLensParameters(self):
        cdef gravitationallens.MultiPlaneContainerParamsPtrConst pParams
        cdef unique_ptr[gravitationallens.GravitationalLens] pSubLens
        cdef double z
        
        self._check()
        pParams = dynamic_cast[gravitationallens.MultiPlaneContainerParamsPtrConst](GravitationalLens._getLens(self).getLensParameters())
        if pParams == NULL:
            raise LensException("Unexpected: parameters are not those of a MultiPlaneContainer")

        params = [ ]
        for i in range(pParams.getNumberOfLenses()):
            z = pParams.getRedshift(i)
            pSubLens = pParams.getLens(i).createCopy()
            if pSubLens.get() == NULL:
                raise LensException("Unexpected: can't create a copy of a contained lens")

            params.append({
                "lens": GravitationalLens._finalizeLoadedLens(pSubLens),
                "z": z
            })
        
        return params

from .privlenses import createLensFromLenstoolFile, createEquivalentPotentialGridLens
