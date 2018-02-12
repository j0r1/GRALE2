
cdef extern from "grale/triangleindices.h" namespace "grale":
    cdef cppclass TriangleIndices:
        int getIndex(int i)

