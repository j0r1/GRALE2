from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from libcpp.memory cimport unique_ptr

cimport grale.serut as serut
cimport grale.errut as errut

cdef extern from "grale/retraceparameters.h" namespace "grale":
    cdef cppclass TraceParameters:
        errut.bool_t read(serut.SerializationInterface &si, unique_ptr[TraceParameters] &params)
        errut.bool_t write(serut.SerializationInterface &si)

cdef extern from "grale/retraceparameters.h" namespace "grale":
    cdef cppclass NoTraceParameters(TraceParameters):
        NoTraceParameters()

cdef extern from "grale/retraceparameters.h" namespace "grale":
    cdef cppclass SingleStepNewtonTraceParams(TraceParameters):
        SingleStepNewtonTraceParams()

ctypedef const SingleStepNewtonTraceParams* SingleStepNewtonTraceParamsPtrConst

cdef extern from "grale/retraceparameters.h" namespace "grale":
    cdef cppclass MultiStepNewtonTraceParams(TraceParameters):
        MultiStepNewtonTraceParams(size_t numEvaluations)
        size_t getNumberOfEvaluations() const

ctypedef const MultiStepNewtonTraceParams* MultiStepNewtonTraceParamsPtrConst

cdef extern from "grale/retraceparameters.h" namespace "grale":
    cdef cppclass ExpandedMultiStepNewtonTraceParams(TraceParameters):
        ExpandedMultiStepNewtonTraceParams(size_t numEvalsPerStartPosition, size_t numMaxGridSteps, double acceptThreshold, double gridSpacing)
        size_t getNumberOfEvaluationsPerStartPosition() const
        size_t getMaximumNumberOfGridSteps() const
        double getAcceptanceThreshold() const
        double getGridSpacing() const

ctypedef const ExpandedMultiStepNewtonTraceParams* ExpandedMultiStepNewtonTraceParamsPtrConst

