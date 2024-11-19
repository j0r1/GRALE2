from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool
from libc.stdint cimport uint64_t

cimport grale.serut as serut
cimport grale.errut as errut
cimport grale.imagesdataextended as imagesdataextended
cimport grale.cppgrid as grid
cimport grale.gravitationallens as gravitationallens
cimport grale.configurationparameters as configurationparameters
cimport grale.vector2d as vector2d
cimport grale.lensinversionbasislensinfo as lensinversionbasislensinfo
cimport grale.scalesearchparameters as scalesearchparameters

cdef extern from "grale/lensinversionparameterssingleplanecpu.h" namespace "grale::LensInversionParametersSinglePlaneCPU":
    cdef enum BasisFunctionType:
        PlummerBasis,
        SquareBasis,
        GaussBasis

    cdef enum MassSheetSearchType:
        NoSheet,
        Genome

cdef extern from "grale/lensinversionparameterssingleplanecpu.h" namespace "grale":
    cdef cppclass LensInversionParametersSinglePlaneCPU(errut.ErrorBase):
        LensInversionParametersSinglePlaneCPU()
        LensInversionParametersSinglePlaneCPU(
                     const vector[shared_ptr[imagesdataextended.ImagesDataExtended]] &images, 
                     const vector[lensinversionbasislensinfo.LensInversionBasisLensInfo] &basisLenses, 
                     double D_d,
                     double z_d,
                     double massscale,
                     bool allowNegativeValues,
                     const gravitationallens.GravitationalLens *pBaseLens,
                     const gravitationallens.GravitationalLens *pSheetLens,
                     const configurationparameters.ConfigurationParameters *pFitnessObjectParams,
                     scalesearchparameters.ScaleSearchParameters &scaleSearchParams,
                     bool randomizeImagePosition,
                     uint64_t initialUncertSeed
                     )

        LensInversionParametersSinglePlaneCPU(
                     const vector[shared_ptr[imagesdataextended.ImagesDataExtended]] &images, 
                     const vector[grid.GridSquare] &gridsquares, 
                     double D_d,
                     double z_d,
                     double massscale,
                     bool useweights,
                     BasisFunctionType basisFunction,
                     bool allowNegativeValues,
                     const gravitationallens.GravitationalLens *pBaseLens,
                     const gravitationallens.GravitationalLens *pSheetLens,
                     const configurationparameters.ConfigurationParameters *pFitnessObjectParams,
                     scalesearchparameters.ScaleSearchParameters &scaleSearchParams,
                     bool randomizeImagePosition,
                     uint64_t initialUncertSeed
                     )

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si)

        @staticmethod
        shared_ptr[gravitationallens.GravitationalLens] createDefaultSheetLens(MassSheetSearchType t, double Dd)

