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

cdef extern from "grale/gridlensinversionparameters.h" namespace "grale::GridLensInversionParameters":
    cdef enum BasisFunctionType:
        PlummerBasis,
        SquareBasis,
        GaussBasis

    cdef enum MassSheetSearchType:
        NoSheet,
        Genome

cdef extern from "grale/gridlensinversionparameters.h" namespace "grale":
    cdef cppclass GridLensInversionParameters(errut.ErrorBase):
        GridLensInversionParameters()
        GridLensInversionParameters(int maxgenerations,
                     const vector[shared_ptr[imagesdataextended.ImagesDataExtended]] &images, 
                     const vector[grid.GridSquare] &gridsquares, 
                     double D_d,
                     double z_d,
                     double massscale,
                     bool useweights,
                     BasisFunctionType basisFunction,
                     bool allowNegativeValues,
                     const gravitationallens.GravitationalLens *pBaseLens,
                     MassSheetSearchType sheetSearchType,
                     const configurationparameters.ConfigurationParameters *pFitnessObjectParams,
                     bool wideSearch)

        bool read(serut.SerializationInterface &si)
        bool write(serut.SerializationInterface &si)

