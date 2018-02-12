from libcpp cimport bool
from libcpp.vector cimport vector

cimport grale.vector2d
ctypedef grale.vector2d.Vector2Dd Vector2Dd

cdef extern from "grale/polygon2d.h" namespace "grale":
    cdef cppclass Polygon2Dd:
        Polygon2Dd()
        void init(vector[Vector2Dd] &points, bool calcHull)
