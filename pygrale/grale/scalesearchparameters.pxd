from libcpp cimport bool

cimport grale.errut as errut

cdef extern from "grale/lensinversionparameterssingleplanecpu.h" namespace "grale":
    cdef cppclass ScaleSearchParameters(errut.ErrorBase):
        ScaleSearchParameters(float startFactor, float stopFactor, int numIt, int firstItSteps, int subseqItSteps)
        ScaleSearchParameters(bool wideSearch)
        ScaleSearchParameters()

        float getStartFactor()
        float getStopFactor()
        int getNumberOfIterations()
        int getStepsOnFirstIteration()
        int getStepsOnSubsequentIterations()
