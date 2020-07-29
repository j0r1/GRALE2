cdef extern from "grale/cosmology.h" namespace "grale":
    cdef cppclass Cosmology:
        Cosmology()
        Cosmology(double h, double Wm, double Wr, double Wv, double w)