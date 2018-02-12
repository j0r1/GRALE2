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

"""

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.cast cimport dynamic_cast
from cython.operator cimport dereference as deref
import cython
import numpy as np
cimport numpy as np
from cpython.array cimport array,clone
import struct
import random
from . import constants
from . import privutil
import copy

cimport grale.gravitationallens as gravitationallens
cimport grale.vector2d as vector2d
cimport grale.serut as serut

include "stringwrappers.pyx"

ctypedef np.ndarray ndarray
ctypedef vector2d.Vector2Dd Vector2Dd

class LensException(Exception):
    """This exception is raised if something goes wrong in the :class:`GravitationalLens` derived classes."""
    pass

# We use a random number as an identifier to make sure that GravitationalLens
# isn't allocated on its own, but only through subclassing
_gravLensRndId = int(random.random()*0xFFFFFFFF)

def getCriticalDensity(double Dd, double Ds, double Dds):
    """getCriticalDensity(Dd, Ds, Dds)
    
    Returns the critical density corresponding to angular diameter distances
    Dd (observer to lens), Ds (observer to source) and Dds (lens to source).

    This is defined as:

    .. code-block:: python

        Sigma_crit = c^2/(4*pi*G*Dd) * (Ds/Ds)

    """
    return constants.SPEED_C**2/(4.0*np.pi*constants.CONST_G*Dd)*(Ds/Dds)

cdef class GravitationalLens:
    """A base class for a gravitational lens."""

    cdef gravitationallens.GravitationalLens *m_pLens

    @staticmethod
    cdef gravitationallens.GravitationalLens * _getLens(GravitationalLens l):
        return l.m_pLens

    def __cinit__(self):
        self.m_pLens = NULL

    def __init__(self, checkId = None):
        if checkId is None or checkId != _gravLensRndId:
            raise LensException("A 'GravitationalLens' object should not be allocated directly, only through subclasses.")

    def __dealloc__(self):
        del self.m_pLens

    cdef _setLens(self, gravitationallens.GravitationalLens *pLens):
        del self.m_pLens
        self.m_pLens = pLens

    cdef _check(self):
        if self.m_pLens == NULL:
            raise LensException("No internal lens object has been set")
    
    cdef _traceTheta1(self, double Ds, double Dds, np.ndarray[double, ndim=1] thetas):
        
        cdef Vector2Dd beta
        cdef np.ndarray[double,ndim=1] betas
        cdef int l,i

        if thetas.shape[0] % 2 == 0:
            l = thetas.shape[0]
            betas = np.zeros([l], dtype = np.double)
            for i in range(0,l,2):
                if not self.m_pLens.traceTheta(Ds, Dds, Vector2Dd(thetas[i], thetas[i+1]), cython.address(beta)):
                    raise LensException(S(self.m_pLens.getErrorString()))
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
                if not self.m_pLens.getAlphaVector(Vector2Dd(thetas[i], thetas[i+1]), cython.address(alpha)):
                    raise LensException(S(self.m_pLens.getErrorString()))
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
                if not self.m_pLens.getAlphaVectorDerivatives(Vector2Dd(thetas[i], thetas[i+1]), axx, ayy, axy):
                    raise LensException(S(self.m_pLens.getErrorString()))
                derivatives[j] = axx
                derivatives[j+1] = ayy
                derivatives[j+2] = axy
                j += 3
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
                invMags[j] = self.m_pLens.getInverseMagnification(Ds, Dds, Vector2Dd(thetas[i], thetas[i+1]))
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
                dens[j] = self.m_pLens.getSurfaceMassDensity(Vector2Dd(thetas[i], thetas[i+1]))
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
                self.m_pLens.getProjectedPotential(Ds, Dds, Vector2Dd(thetas[i], thetas[i+1]), cython.address(value))
                potential[j] = value
                j += 1
        else:
            raise LensException("Bad 1D array dimensions, must be a multiple of two")

        return potential

    cdef _getRadialMassProfile1(self, np.ndarray[double, ndim=1] thetaRadii):

        cdef int i,l
        cdef double value = 0
        cdef np.ndarray[double,ndim=1] profile
        cdef gravitationallens.SymmetricLens *pSymLens = gravitationallens.SymmetricLens.cast(self.m_pLens)
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
        cdef gravitationallens.SymmetricLens *pSymLens = gravitationallens.SymmetricLens.cast(self.m_pLens)
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
        return self._reshapeAndCall1D(lambda x : self._traceTheta1(float(Ds), float(Dds), x), thetas, 2, 2)

    def getAlphaVector(self, thetas):
        """getAlphaVector(thetas)
        
        Returns the deflection angles at the positions in 'thetas'.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getAlphaVector1(x), thetas, 2, 2)

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
        return self._reshapeAndCall1D(lambda x : self._getAlphaVectorDerivatives1(x), thetas, 2, 3)

    def getInverseMagnification(self, Ds, Dds, thetas):
        """getInverseMagnification(Ds, Dds, thetas)
        
        Returns the inverse magnification values at the positions in 'thetas', for a source
        plane with angular diameter distances Ds and Dds.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getInverseMagnification1(float(Ds), float(Dds), x), thetas, 2, 1)

    def getSurfaceMassDensity(self, thetas):
        """getSurfaceMassDensity(thetas)
        
        Returns the surface mass density of the gravitational lens at the positions
        in 'thetas'.
        """
        thetas = np.array(thetas)
        return self._reshapeAndCall1D(lambda x : self._getSurfaceMassDensity1(x), thetas, 2, 1)

    def getSurfaceMassDensityMap(self, bottomLeft, topRight, int numX, int numY, feedbackObject = "default", renderer = None, reduceToPixels = True):
        """getSurfaceMassDensityMap(bottomLeft, topRight, numX, numY, feedbackObject = "default", renderer = None, reduceToPixels = True)
        
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
                    massMap[yi,xi] = self.m_pLens.getSurfaceMassDensity(Vector2Dd(x,y))

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
        return self._reshapeAndCall1D(lambda x : self._getProjectedPotential1(float(Ds), float(Dds), x), thetas, 2, 1)

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
        return self.m_pLens.getLensDistance()

    def setLensDistance(self, double Dd):
        """setLensDistance(Dd)
        
        This function changes the angular diameter distance to the lens to Dd.
        """
        self._check()
        if Dd <= 0:
            raise LensException("Invalid distance")
        self.m_pLens.setLensDistance(Dd)

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
        factor = self.m_pLens.getLensDistance() * (1.0 + zd) / constants.SPEED_C * (Ds/Dds)
        delay = (0.5*distances - potential)*factor
        return delay

    def setDerivativeAngularDistanceScale(self, double scale):
        """setDerivativeAngularDistanceScale(scale)
        
        In case the :func:`getAlphaVectorDerivatives` function is not provided
        for a specific lens, the derivatives will be calculated numerically
        using points that are this distance apart.
        """
        self._check()
        self.m_pLens.setDerivativeAngularDistanceScale(scale)

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
        cdef gravitationallens.GravitationalLens *pLens = NULL
        if not gravitationallens.GravitationalLens.load(B(fileName), cython.address(pLens), errorString):
            raise LensException(S(errorString))

        return GravitationalLens._finalizeLoadedLens(pLens)

    @staticmethod
    cdef _finalizeLoadedLens(gravitationallens.GravitationalLens *pLens):
        cdef gravitationallens.LensType t = pLens.getLensType()
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
        else: # Unknown, can still use the interface
            l = GravitationalLens(_gravLensRndId)

        l.m_pLens = pLens
        return l

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(bytes b)
        
        Instantiates a gravitational lens object from the byte array previously generated
        by :func:`toBytes`.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef serut.MemorySerializer *m = new serut.MemorySerializer(buf.data.as_voidptr, len(b), NULL, 0)
        cdef gravitationallens.GravitationalLens *pLens = NULL
        cdef string errorString
        
        try:
            if not gravitationallens.GravitationalLens.read(deref(m), cython.address(pLens), errorString):
                raise LensException(S(errorString))

            return GravitationalLens._finalizeLoadedLens(pLens)
        finally:
            del m

    def save(self, fileName):
        """save(fileName)
        
        Write this gravitational lens object to the file with name 'fileName'.
        Can be loaded again using :func:`load`.
        """
        self._check()
        if not self.m_pLens.save(B(fileName)):
            raise LensException(S(self.m_pLens.getErrorString()))

    def toBytes(self):
        """toBytes()
        
        Converts the lens object to a byte array, which can be converted to a real
        lens object again using :func:`fromBytes`
        """
        cdef serut.VectorSerializer vSer

        self._check()
        if not self.m_pLens.write(vSer):
            raise LensException(S(self.m_pLens.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

    # Will be overridden
    cdef gravitationallens.GravitationalLens* _allocLens(self) except NULL:
        return NULL

    # Will be overridden
    cdef gravitationallens.GravitationalLensParams *_allocParams(self, params) except NULL:
        return NULL

    def _lensInit(self, Dd, params):
        cdef gravitationallens.GravitationalLensParams *lp = NULL
        cdef gravitationallens.GravitationalLens *pLens = NULL

        if Dd is None and params is None: # Don't initialize yet
            return

        lp = self._allocParams(params)
        pLens = self._allocLens()
        
        if not pLens.init(Dd, lp):
            err = pLens.getErrorString()
            del pLens
            del lp
            raise LensException(S(err))
        
        del lp
        self._setLens(pLens)

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
        cdef gravitationallens.GravitationalLens *pSubLens
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
            if pSubLens == NULL:
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


