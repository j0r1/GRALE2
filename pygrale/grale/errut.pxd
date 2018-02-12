from libcpp.string cimport string

cdef extern from "errut/errorbase.h" namespace "errut":

    cdef cppclass ErrorBase:
        string getErrorString()

