from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from libcpp.memory cimport unique_ptr

cimport grale.serut as serut
cimport grale.errut as errut

cdef extern from "grale/eaparameters.h" namespace "grale":
    cdef cppclass EAParameters:
        errut.bool_t read(serut.SerializationInterface &si, unique_ptr[EAParameters] &params)
        errut.bool_t write(serut.SerializationInterface &si)

cdef extern from "grale/gaparameters.h" namespace "grale":
    cdef cppclass GAParameters(EAParameters):
        GAParameters(double selectionPressure, cbool useElitism, cbool alwaysIncludeBest, double crossoverRate)

        double getSelectionPressure() const
        cbool getUseElitism() const
        cbool getAlwaysIncludeBest() const
        double getCrossOverRate() const

ctypedef const GAParameters* GAParametersPtrConst

cdef extern from "grale/deparameters.h" namespace "grale":
    cdef cppclass JADEParameters(EAParameters):
        JADEParameters()
