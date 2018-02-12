from libcpp cimport bool as cbool
from libcpp.vector cimport vector
cimport grale.vector2d as vector2d

cdef extern from "grale/contourfinder.h" namespace "grale":
    cdef cppclass ContourFinder:
        ContourFinder(vector[double] &values, vector2d.Vector2Dd bl, vector2d.Vector2Dd tr, int numX, int numY)

        vector[vector[vector2d.Vector2Dd]] findContour(double level)
