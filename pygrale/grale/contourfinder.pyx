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
    pass

cdef class ContourFinder:

    cdef cppcontourfinder.ContourFinder *m_pContourFinder

    def __cinit__(self):
        self.m_pContourFinder = NULL

    def __dealloc__(self):
        del self.m_pContourFinder

    def __init__(self, np.ndarray[double, ndim=2] values, bottomLeft, topRight):
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
        
        self.m_pContourFinder = new cppcontourfinder.ContourFinder(valueVector, bl, tr, numX, numY)

    def _check(self):
        if self.m_pContourFinder == NULL:
            raise ContourFinderException("No internal contour finder instance was allocated")

    def findContour(self, level):
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
