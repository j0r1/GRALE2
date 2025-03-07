from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool
from libc.stdint cimport uint64_t
from libcpp.pair cimport pair

cimport grale.serut as serut
cimport grale.errut as errut
cimport grale.imagesdataextended as imagesdataextended
cimport grale.cppgrid as grid
cimport grale.gravitationallens as gravitationallens
cimport grale.configurationparameters as configurationparameters
cimport grale.vector2d as vector2d
cimport grale.lensinversionbasislensinfo as lensinversionbasislensinfo
cimport grale.scalesearchparameters as scalesearchparameters

cdef extern from "grale/lensinversionparametersparametricsingleplane.h" namespace "grale":
    cdef cppclass LensInversionParametersParametricSinglePlane(errut.ErrorBase):
        LensInversionParametersParametricSinglePlane()
        LensInversionParametersParametricSinglePlane(
                const vector[shared_ptr[imagesdataextended.ImagesDataExtended]] &images,
                double Dd, double zd,
                const gravitationallens.GravitationalLens &templateLens,
                double deflectionScale, double potentialScale,
                const vector[int] &offsets,
                const vector[float] &initMin, const vector[float] &initMax,
                const vector[float] &hardMin, const vector[float] &hardMax,
                bool infOnBoundsViolation,
                const configurationparameters.ConfigurationParameters &fitnessObjectParams,
                int devIdx,
                bool randomizeImagePositions,
                uint64_t initialUncertSeed,
                const vector[pair[size_t,string]] &originParameterMapping,
                size_t numOriginParams,
                bool allowUnusedPriors,
                const vector[bool] &retraceImages,
                size_t numRetraceSteps,
                double sourcePlaneDistThreshold,
                string clPriorCode,
                bool allowEqualInitRange
                )

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si)

