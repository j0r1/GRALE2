from libcpp cimport bool as cbool
from libcpp.vector cimport vector
from libcpp.string cimport string
cimport grale.vector2d as vector2d

cdef extern from "grale/contourfinder.h" namespace "grale":
    cdef cppclass ContourFinder:
        ContourFinder(vector[double] &values, vector2d.Vector2Dd bl, vector2d.Vector2Dd tr, int numX, int numY)

        vector[vector[vector2d.Vector2Dd]] findContour(double level)

cdef extern from "multicontourfinder.h":
    cdef cppclass MultiContourFinder(ContourFinder):
        MultiContourFinder(vector[double] &values, vector2d.Vector2Dd bl, vector2d.Vector2Dd tr, int numX, int numY)

        string getErrorString()
        cbool findContours(vector[double] &level, int numThreads)
        vector[vector[vector[vector2d.Vector2Dd]]] &getContours()

