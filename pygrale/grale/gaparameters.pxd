from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool

cimport grale.serut as serut
cimport grale.errut as errut

cdef extern from "grale/gaparameters.h" namespace "grale":
    cdef cppclass GAParameters(errut.ErrorBase):
        GAParameters(double selectionPressure, cbool useElitism, cbool alwaysIncludeBest, double crossoverRate)

        cbool read(serut.SerializationInterface &si)
        cbool write(serut.SerializationInterface &si)

        double getSelectionPressure() const
        cbool getUseElitism() const
        cbool getAlwaysIncludeBest() const
        double getCrossOverRate() const


