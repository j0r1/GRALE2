"""This module defines a :class:`ContourFinder` class, with which
you can look for specific contour levels in a NumPy grid. This code
is used to look for critical lines for example, as the contour where
the inverse magnification is zero. It can also be used to search for
specific contour levels in the mass map of a gravitational lens."""

from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cython.operator cimport dereference as deref
import cython
import numpy as np
cimport numpy as np

include "stringwrappers.pyx"

ctypedef np.ndarray ndarray

cimport grale.cppcontourfinder as cppcontourfinder
from grale.vector2d cimport Vector2Dd

class ContourFinderException(Exception):
    """An exception that will be raised if something goes wrong when
    looking for a contour."""
    pass

cdef class ContourFinder:
    """After initialization with a specific grid, you can use the
    :func:`findContour` member function to calculate a specific
    contour."""

    cdef cppcontourfinder.MultiContourFinder *m_pContourFinder

    def __cinit__(self):
        self.m_pContourFinder = NULL

    def __dealloc__(self):
        del self.m_pContourFinder

    def __init__(self, np.ndarray[double, ndim=2] values, bottomLeft, topRight):
        """__init__(values, bottomleft, topright)

        Initializes the instance so that a specific grid is used to look
        for contours.

        Arguments:
        
         - `values`: a 2D NumPy grid of values, which will be interpreted as a
           height map.
         - `bottomLeft`: the value at [0,0] in `values` will correspond to this
           (X,Y) coordinate.
         - `topRight`: if the shape of `values` is (Ny, Nx), then this is the
           (X,Y) coordinate of the pixel [Ny-1, Nx-1]
        """
        cdef int x, y, numX, numY, offset
        cdef Vector2Dd bl = Vector2Dd(bottomLeft[0], bottomLeft[1])
        cdef Vector2Dd tr = Vector2Dd(topRight[0], topRight[1])
        cdef vector[double] valueVector

        if values is None:
            raise ContourFinderException("Not values were specified")

        numY, numX = values.shape[0], values.shape[1]

        if numY < 2 or numX < 2:
            raise ContourFinderException("The grid must consist of at least two points in each direction")

        valueVector.resize(numX*numY)
        for y in range(numY):
            offset = y*numX
            for x in range(numX):
                valueVector[x+offset] = values[y,x]
        
        self.m_pContourFinder = new cppcontourfinder.MultiContourFinder(valueVector, bl, tr, numX, numY)

    def _check(self):
        if self.m_pContourFinder == NULL:
            raise ContourFinderException("No internal contour finder instance was allocated")

    def findContour(self, level):
        """findContour(level)
        
        Within the height map as specified during initialization, look for the contours
        that correspond to the height `level`. As several distinct lines may be found, the
        function returns a list of lists, where the second list contains the consecutive coordinates of 
        a part of the contour."""
        cdef vector[vector[Vector2Dd]] contour
        cdef int i, j

        self._check()
        contour = self.m_pContourFinder.findContour(level)

        contourParts = [ ]
        for i in range(int(contour.size())):
            part = [ ]
            for j in range(int(contour[i].size())):
                part.append([contour[i][j].getX(), contour[i][j].getY()])

            contourParts.append(np.array(part))

        return contourParts

    def findMultipleContours(self, levels, numThreads = -1):
        """findMultipleContours(levels, numThreads = -1)
        Similar to :func:`findContour`, but for several levels. If `numThreads` is positive,
        this is the amount of contours that will be calculated in parallel, speeding things
        up if you have multiple cores. If set to 0 or a negative value, the number of threads
        will be set to the value detected by the `multiprocessing.cpucount <https://docs.python.org/3/library/multiprocessing.html>`_
        function."""
        cdef vector[double] vecLevels
        cdef vector[vector[vector[Vector2Dd]]] contours

        self._check()

        if numThreads <= 0:
            import multiprocessing
            numThreads = multiprocessing.cpu_count()

        for l in levels:
            vecLevels.push_back(l)

        if not self.m_pContourFinder.findContours(vecLevels, numThreads):
            raise ContourFinderException(S(self.m_pContourFinder.getErrorString()))

        contours = self.m_pContourFinder.getContours()

        levelData = [ ]
        for l in range(int(contours.size())):
            contourParts = [ ]
            for i in range(int(contours[l].size())):
                part = [ ]
                for j in range(int(contours[l][i].size())):
                    part.append([contours[l][i][j].getX(), contours[l][i][j].getY()])

                contourParts.append(np.array(part))

            levelData.append(contourParts)

        return levelData
