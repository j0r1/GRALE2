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
        GAParameters(double selectionPressure, cbool useElitism, cbool alwaysIncludeBest, double crossoverRate, double smallMutationSize)

        double getSelectionPressure() const
        cbool getUseElitism() const
        cbool getAlwaysIncludeBest() const
        double getCrossOverRate() const
        double getSmallMutationSize() const

ctypedef const GAParameters* GAParametersPtrConst

cdef extern from "grale/deparameters.h" namespace "grale":
    cdef cppclass DEParameters(EAParameters):
        DEParameters(double F, double CR, cbool needStrictlyBetter)
        double getF() const
        double getCR() const
        cbool getNeedStrictlyBetter() const

ctypedef const DEParameters* DEParametersPtrConst

cdef extern from "grale/deparameters.h" namespace "grale":
    cdef cppclass JADEParameters(EAParameters):
        JADEParameters(double p, double c, cbool useArchive, double initMuF, double initMuCR, cbool needStrictlyBetter)
        double getBestFraction_p() const
        double getParameterUpdateFraction_c() const
        cbool useExternalArchive() const
        double getInitialMeanF() const
        double getInitialMeanCR() const
        cbool getNeedStrictlyBetter() const

ctypedef const JADEParameters* JADEParametersPtrConst

