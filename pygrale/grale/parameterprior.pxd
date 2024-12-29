from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.map cimport map as cmap
from libcpp.memory cimport shared_ptr, unique_ptr
from libcpp cimport bool

cimport grale.serut as serut
cimport grale.errut as errut

cdef extern from "grale/parameterprior.h" namespace "grale":
    cdef cppclass ParameterPrior:
        pass

cdef extern from "grale/parameterprior.h" namespace "grale":
    cdef cppclass NoParameterPrior(ParameterPrior):
        NoParameterPrior()

cdef extern from "grale/parameterprior.h" namespace "grale":
    cdef cppclass GaussianParameterPrior(ParameterPrior):
        GaussianParameterPrior(float mu, float sigma)
