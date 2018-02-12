from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
import cython
import numpy as np
cimport numpy as np

include "stringwrappers.pyx"

ctypedef np.ndarray ndarray

cimport grale.cppgridfunction as cppgridfunction
from grale.vector2d cimport Vector2Dd

class ResampleException(Exception):
    pass

# TODO: can we replace this using a gridfunction?
def resample2DArray(np.ndarray[double, ndim=2] origMap, int newDim1, int newDim2):
    cdef np.ndarray[double,ndim=2] resampMap
    cdef double value
    cdef int x1, x2, y1, y2
    cdef int y1fraci, y2fraci
    cdef double y1frac, y2frac
    cdef double y1d, y2d
    cdef double v11, v12, v21, v22

    if origMap is None:
        raise ResampleException("No source map was specified")
    if origMap.shape[0] < 2 or origMap.shape[1] < 2:
        raise ResampleException("The source map should have at least two points in every direction")
    if newDim1 < 2 or newDim2 < 2:
        raise ResampleException("Both target map dimensions should be at least two")

    resampMap = np.zeros([newDim1, newDim2], dtype = np.double)
    for x1 in range(0, newDim1):
        for x2 in range(0, newDim2):

            if x1 == 0:
                y1 = 0
                y1fraci = 0
                y1frac = 0
            elif x1 == newDim1-1:
                y1 = origMap.shape[0]-1
                y1fraci = 0
                y1frac = 0
            else:
                y1fraci = 1
                y1d = (<double>x1/<double>(newDim1-1))*<double>(origMap.shape[0]-1)
                y1 = <int>y1d
                y1frac = y1d - <double>y1

            if x2 == 0:
                y2 = 0
                y2fraci = 0
                y2frac = 0
            elif x2 == newDim2-1:
                y2 = origMap.shape[1]-1
                y2fraci = 0
                y2frac = 0
            else:
                y2fraci = 1
                y2d = (<double>x2/<double>(newDim2-1))*<double>(origMap.shape[1]-1)
                y2 = <int>y2d
                y2frac = y2d - <double>y2
            
            v11 = origMap[y1, y2]

            if y1fraci:
                if y2fraci:
                    v22 = origMap[y1+1, y2+1]
                    v12 = origMap[y1, y2+1]
                    v21 = origMap[y1+1, y2]
                else:
                    v22 = origMap[y1+1, y2]
                    v12 = origMap[y1,   y2]
                    v21 = origMap[y1+1, y2]
            else:
                if y2fraci:
                    v22 = origMap[y1, y2+1]
                    v12 = origMap[y1, y2+1]
                    v21 = origMap[y1, y2]
                else:
                    v22 = origMap[y1, y2]
                    v12 = origMap[y1, y2]
                    v21 = origMap[y1, y2]

            #print "v11 = %g v12 = %g v21 = %g v22 = %g y1frac = %g y2frac = %g" % (v11,v12,v21,v22,y1frac,y2frac)
            value = v11*(1.0-y1frac)*(1.0-y2frac) + v21*y1frac*(1.0-y2frac) + v12*(1.0-y1frac)*y2frac + v22*y1frac*y2frac
            resampMap[x1,x2] = value

    return resampMap

