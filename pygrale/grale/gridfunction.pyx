"""This module defines a :class:`GridFunction` with which you can
treat a 2D NumPy array as a function"""
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from libc.math cimport sin, cos
from cython.operator cimport dereference as deref
import cython
import math
import numpy as np
cimport numpy as np

include "stringwrappers.pyx"

ctypedef np.ndarray ndarray

cimport grale.cppgridfunction as cppgridfunction
from grale.vector2d cimport Vector2Dd

class GridFunctionException(Exception):
    """An exception that will be thrown if something goes wrong in
    the :class:`GridFunction` class."""
    pass

cdef class GridFunction:
    """An instance of this class can be used to calculate values at specific
    points using the :func:`evaluate` member function. The calculated values
    are based on the values of a NumPy grid, which are said to correspond to
    a specific coordinate region.

    The constructor itself allows this grid to be rotated as well, but the
    member functions :func:`createFromCorners` and :func:`createFromFITS`
    may be easier to use if the input values are not aligned with the
    coordinate axes.
    """

    cdef cppgridfunction.GridFunction *m_pGridFunction
    cdef vector[double] m_values
    cdef double m_rotAngleRadians

    def __cinit__(self):
        self.m_pGridFunction = NULL
        self.m_rotAngleRadians = 0.0

    # Region is rotated counterclockwise over rotationAngleRadians around bottomLeft
    def __init__(self, np.ndarray[double, ndim=2] values, bottomLeft, topRight, cbool asPixels = False, double rotationAngleRadians = 0.0):
        """__init__(values, bottomLeft, topRight, asPixels = False, rotationAngleRadians = 0.0)

        Constructs an instance to correspond with a specific grid of values.

        Arguments:

         - `values`: a 2D NumPy array describing the values this function should be based on.
         - `bottomLeft` and `topRight`: describe the area in the coordinate space that the
           grid values are mapped to.
         - `asPixels`: if ``False``, the [0,0] and [Ny-1, Nx-1] values of the grid, describe
           the values precisely at the `bottomLeft` and `topRight` X,Y coordinates and each
           grid value will correspond to a specific point. If a different coordinate is used
           int the :func:`evaluate` function, then bi-linear interpolation is used.
           
           In case `asPixels` is ``True``, then the area is divided into pixels which have a
           constant value inside, based on the grid values of course.
         - `rotationAngleRadians`: if specified, the region that's obtained from the
           `bottomLeft` and `topRight` X,Y coordinates will be rotated further by this angle.
           Note however that it may be easier to use the :func:`createFromCorners` static
           function if you need to specify an area that's not aligned with the coordinate
           axes.
        """

        cdef int x, y, numY, numX, offset
        cdef Vector2Dd bl = Vector2Dd(bottomLeft[0], bottomLeft[1])
        cdef Vector2Dd tr = Vector2Dd(topRight[0], topRight[1])

        if values is None:
            raise GridFunctionException("No values were specified")

        numY, numX = values.shape[0], values.shape[1]

        if numY < 2 or numX < 2:
            raise GridFunctionException("The grid must consist of at least two points in each direction")

        self.m_values.resize(numX*numY)
        for y in range(numY):
            offset = y*numX
            for x in range(numX):
                self.m_values[x+offset] = values[y,x]
        
        self.m_pGridFunction = new cppgridfunction.GridFunction(<double*>&(self.m_values[0]), bl, tr, numX, numY, asPixels)
        self.m_rotAngleRadians = rotationAngleRadians

    @staticmethod
    def createFromCorners(values, p0, p1, p2, asPixels = False):
        """createFromCorners(values, p0, p1, p2, asPixels = False)
        
        Similar to the contructor, an instance of the grid function is based on
        `values`, but instead of `bottomLeft` and `topRight` corners, this function
        requires you to specify the corner `p0` that corresponds with the (0,0) pixel,
        the corner `p1` that corresponds to the (0, Nx-1) corner and `p2` that
        corresponds to (Ny-1, 0).
        """
        
        getDiff = lambda a, b: (b[0]-a[0],b[1]-a[1])
        dotProd = lambda a, b: a[0]*b[0] + a[1]*b[1]
        getLength = lambda a: dotProd(a,a)**0.5
        def normalize(a):
            l = getLength(a)
            return (a[0]/l, a[1]/l)

        d1 = getDiff(p0, p1)
        d2 = getDiff(p0, p2)

        if abs(dotProd(normalize(d1), normalize(d2))) > 1e-6:
            raise GridFunctionException("Specified grid vectors are not rectangular")

        # Find the angle that d1 makes
        n1 = normalize(d1)
        angle = math.acos(n1[0])
        if n1[1] < 0:
            angle = -angle

        # Rotate d1 and d2 clockwise over angle so that they're aligned again
        def rotate(v, a):
            return (v[0]*math.cos(a) + v[1]*math.sin(a),
                   -v[0]*math.sin(a) + v[1]*math.cos(a))

        D1 = rotate(d1, angle)
        D2 = rotate(d2, angle)

        bl = p0
        tr = [ p0[0] + D1[0], p0[1] + D2[1] ]

        return GridFunction(values, bl, tr, asPixels = asPixels, rotationAngleRadians = angle)
        
    @staticmethod
    def createFromFITS(fitsHDUItem, centerRaDec, asPixels = False):
        """createFromFITS(fitsHDUItem, centerRaDec, asPixels = False)
        
        If you want to use a grid function based on the values in a FITS file,
        you can use this function.

        Arguments:

         - `fitsHDUItem`: the part of the FITS file that you want to use
         - `centerRaDec`: re-center the coordinates so that this RA,Dec value
           becomes the (0,0) coordinate.
         - `asPixels`: same as in the constructor, will determine whether or
           not values are interpolated or constant within a pixel.
        """

        from .images import centerOnPosition
        from .constants import ANGLE_DEGREE
        from astropy import wcs
        
        w = wcs.WCS(fitsHDUItem.header)
        data = fitsHDUItem.data.astype(np.double)
        p0 = np.array(w.all_pix2world(0.5, 0.5, 0))*ANGLE_DEGREE
        p1 = np.array(w.all_pix2world(data.shape[1]-0.5, 0.5, 0))*ANGLE_DEGREE
        p2 = np.array(w.all_pix2world(0.5, data.shape[0]-0.5, 0))*ANGLE_DEGREE
        
        p0 = centerOnPosition(p0, centerRaDec)
        p1 = centerOnPosition(p1, centerRaDec)
        p2 = centerOnPosition(p2, centerRaDec)
        
        return GridFunction.createFromCorners(data, p0, p1, p2, asPixels)

    def __dealloc__(self):
        del self.m_pGridFunction

    def _check(self):
        if self.m_pGridFunction == NULL:
            raise GridFunctionException("No internal grid function has been set")

    cdef _evaluate1D(self, np.ndarray[double, ndim=1] points, cbool _checkRegion):

        cdef int i, j, l
        cdef double value
        cdef np.ndarray[double,ndim=1] values
        cdef Vector2Dd point, diff
        cdef cppgridfunction.GridFunction *pFunc = self.m_pGridFunction
        cdef cbool checkRegion = _checkRegion
        cdef Vector2Dd bl = pFunc.getBottomLeft()
        cdef Vector2Dd tr = pFunc.getTopRight()
        cdef double x0 = bl.getX()
        cdef double x1 = tr.getX()
        cdef double y0 = bl.getY()
        cdef double y1 = tr.getY()
        cdef cbool shouldRotate = True if self.m_rotAngleRadians != 0.0 else False
        cdef double angle = self.m_rotAngleRadians

        if x0 > x1:
            x0, x1 = x1, x0
        if y0 > y1:
            y0, y1 = y1, y0

        if points.shape[0]%2 != 0:
            raise GridFunctionException("Bad 1D array dimensions, must be a multiple of two")

        j = 0
        l = points.shape[0]
        values = np.zeros([l//2], dtype=np.double)
        for i in range(0, l, 2):
            point = Vector2Dd(points[i], points[i+1])
            if shouldRotate:
                diff = point
                diff -= bl
                point = bl
                point +=  Vector2Dd( diff.getX()*cos(angle) + diff.getY()*sin(angle),
                                     -diff.getX()*sin(angle) + diff.getY()*cos(angle))

            if checkRegion:
                if point.getX() < x0 or point.getX() > x1 or point.getY() < y0 or point.getY() > y1:
                    raise GridFunctionException("Not all coordinates are in the region of the GridFunction (failed on ({},{}))".format(point.getX(), point.getY()))
            values[j] = (deref(pFunc))(point)
            j += 1

        return values

    def evaluate(self, points, checkRegion = True):
        """evaluate(points, checkRegion = True)
        
        For all points in `points`, calculate the value based on the specified
        grid. If `checkRegion` is ``True``, an exception will be generated if a point
        lies outside the grid region, otherwise the value will be based on the border
        values.
        """

        cdef int i, l, totalElements

        if points is None:
            raise GridFunctionException("No points were specified")

        self._check()
        
        l = len(points.shape)
        if l == 1:
            return self._evaluate1D(points, checkRegion)

        if l > 1 and points.shape[l-1] == 2:
            totalElements = 1
            for i in range(l):
                totalElements *= points.shape[i]

            outShape = points.shape[:-1] 
            return np.reshape(self._evaluate1D(np.reshape(points,[totalElements]), checkRegion), outShape)

        raise GridFunctionException("Bad points array dimensions")



