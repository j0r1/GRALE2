from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool

cimport grale.serut as serut
cimport grale.errut as errut
cimport grale.imagesdataextended as imagesdataextended
cimport grale.cppgrid as grid
cimport grale.gravitationallens as gravitationallens
cimport grale.configurationparameters as configurationparameters
cimport grale.vector2d as vector2d
cimport grale.lensinversionbasislensinfo as lensinversionbasislensinfo
cimport grale.scalesearchparameters as scalesearchparameters
cimport grale.cppcosmology as cppcosmology

cdef extern from "grale/lensinversionparametersmultiplanegpu.h" namespace "grale":
    cdef cppclass LensInversionParametersMultiPlaneGPU(errut.ErrorBase):
        LensInversionParametersMultiPlaneGPU()
        LensInversionParametersMultiPlaneGPU(const cppcosmology.Cosmology &cosmology,
            const vector[double] &lensRedshifts,
            const vector[vector[shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo]]] &basisLenses,
            const vector[shared_ptr[gravitationallens.GravitationalLens]] &baseLensesPerPlane,
            const vector[shared_ptr[imagesdataextended.ImagesDataExtended]] &sourceImages,
            double massEstimate,
            bool useMassSheets,
            const configurationparameters.ConfigurationParameters *pFitnessObjectParams,
            bool allowNegativeWeights,
            const scalesearchparameters.ScaleSearchParameters &scaleSearchParams,
            int deviceIndex)

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si)
