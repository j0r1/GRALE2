from libcpp.string cimport string
from libcpp cimport bool as cbool

cdef extern from "errut/errorbase.h" namespace "errut":

    cdef cppclass ErrorBase:
        string getErrorString()

cdef extern from "errut/booltype.h" namespace "errut":
    cdef cppclass bool_t:
        string getErrorString()
        cbool success()
