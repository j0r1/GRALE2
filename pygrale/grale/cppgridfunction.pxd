from libcpp cimport bool as cbool
cimport grale.vector2d as vector2d

cdef extern from "grale/gridfunction.h" namespace "grale":
    cdef cppclass GridFunction:
        GridFunction(const double *pValues, vector2d.Vector2Dd bl, vector2d.Vector2Dd tr, int numX, int numY, cbool asPixels)
        double operator()(vector2d.Vector2Dd pos)

        vector2d.Vector2Dd getBottomLeft()
        vector2d.Vector2Dd getTopRight()
