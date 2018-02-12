
cdef extern from "grale/vector2d.h" namespace "grale":
    cdef cppclass Vector2Dd:
        Vector2Dd()
        Vector2Dd(double x, double y)
        double getX()
        double getY()

