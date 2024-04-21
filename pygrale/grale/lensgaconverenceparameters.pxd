from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport grale.serut as serut

cdef extern from "grale/lensgaconvergenceparameters.h" namespace "grale":
    cdef cppclass LensGAConvergenceParameters:
        LensGAConvergenceParameters()

        void setMaximumNumberOfGenerations(size_t s)
        void setHistorySize(size_t s)
        void setConvergenceFactor(double convFactor)

        size_t getMaximumNumberOfGenerations() const
        size_t getConvergenceHistorySize() const
        double getConvergenceFactor() const

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si) const

        string getErrorString()
