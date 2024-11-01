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

cdef extern from "grale/eaparameters.h" namespace "grale":
    cdef cppclass EATestParameters(EAParameters):
        EATestParameters()

ctypedef const EATestParameters* EATestParametersPtrConst

cdef extern from "grale/gaparameters.h" namespace "grale":
    cdef cppclass GAParameters(EAParameters):
        GAParameters(double selectionPressure, cbool useElitism, cbool alwaysIncludeBest, double crossoverRate, double smallMutationSize)

        double getSelectionPressure() const
        cbool getUseElitism() const
        cbool getAlwaysIncludeBest() const
        double getCrossOverRate() const
        double getSmallMutationSize() const

ctypedef const GAParameters* GAParametersPtrConst

cdef extern from "grale/nsga2parameters.h" namespace "grale":
    cdef cppclass NSGA2Parameters(EAParameters):
        NSGA2Parameters(double smallMutationSize)
        double getSmallMutationSize() const

ctypedef const NSGA2Parameters* NSGA2ParametersPtrConst

cdef extern from "grale/nsga2parameters.h" namespace "grale":
    cdef cppclass NSGA2DELikeCrossoverParameters(EAParameters):
        NSGA2DELikeCrossoverParameters(cbool extraParent, float F, float CR)
        cbool useExtraParent()
        float getF()
        float getCR()

ctypedef const NSGA2DELikeCrossoverParameters* NSGA2DELikeCrossoverParametersPtrConst

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

cdef extern from "grale/rndparameters.h" namespace "grale":
    cdef cppclass RNDParameters(EAParameters):
        RNDParameters(double scale)
        double getScale() const

ctypedef const RNDParameters* RNDParametersPtrConst

cdef extern from "grale/mcmcparameters.h" namespace "grale":
    cdef cppclass MCMCParameters(EAParameters):
        MCMCParameters(double a, const string &sampleFileName, size_t sampleGenerations,
                       size_t burnInGenerations,
                       size_t annealGenerationsScale, double alpha0, double alphaMax)

        double getGoodmanWeare_a() const
        const string getSamplesFilename()
        size_t getSampleGenerations() const
        size_t getBurnInGenerations() const
        size_t getAnnealGenerationsTimeScale()
        double getAnnealAlpha0() const
        double getAnnealAlphaMax() const

ctypedef const MCMCParameters* MCMCParametersPtrConst

