cimport grale.vector2d
ctypedef grale.vector2d.Vector2Dd Vector2Dd

cdef extern from "grale/grid.h" namespace "grale":
    cdef cppclass GridSquare:
        GridSquare(Vector2Dd center, double size)
