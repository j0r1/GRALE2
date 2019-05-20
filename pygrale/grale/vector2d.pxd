
cdef extern from "grale/vector2d.h" namespace "grale":
    cdef cppclass Vector2Dd:
        Vector2Dd()
        Vector2Dd(double x, double y)
        double getX()
        double getY()

    cdef cppclass Vector2Df:
        Vector2Df()
        Vector2Df(double x, double y)
        double getX()
        double getY()


