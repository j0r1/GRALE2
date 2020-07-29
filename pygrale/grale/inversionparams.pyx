"""This module contains classes which will contain various settings
for the lens inversion method. They form the bridge between the Python
side of the code, and the C++ one.
"""
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr, make_shared
from libcpp cimport bool as cbool
from cpython.array cimport array,clone
import cython
from cython.operator cimport dereference as deref
from . import images
from . import lenses

cimport grale.configurationparameters as configurationparameters
cimport grale.lensinversionparameterssingleplanecpu as lensinversionparameterssingleplanecpu
cimport grale.lensinversionparametersmultiplanegpu as lensinversionparametersmultiplanegpu
cimport grale.serut as serut
cimport grale.imagesdataextended as imagesdataextended
cimport grale.cppgrid as grid
cimport grale.vector2d as vector2d
cimport grale.gravitationallens as gravitationallens
cimport grale.configurationparameters as configurationparameters
cimport grale.gaparameters as gaparameters
cimport grale.lensinversionbasislensinfo as lensinversionbasislensinfo
cimport grale.scalesearchparameters as scalesearchparameters
cimport grale.cppcosmology as cppcosmology

include "stringwrappers.pyx"

class ConfigurationParametersException(Exception):
    """This exception is raised in case somethings goes wrong in the
    :class:`ConfigurationParameters` class."""
    pass

cdef class ConfigurationParameters(object):
    """This class describes configuration parameters for a specific lens inversion
    module. It is the Python link to the similarly named C++ class."""

    cdef configurationparameters.ConfigurationParameters *m_pParams

    def __cinit__(self):
        self.m_pParams = new configurationparameters.ConfigurationParameters()

    @staticmethod
    cdef configurationparameters.ConfigurationParameters * _getConfigurationParameters(confParams):
        return confParams.m_pParams

    def __init__(self, parameterDict = None):
        """__init__(parameterDict = None)

        Initializes this object, storing the key/value entries from `parameterDict`
        (if present). Note that not all dictionaries are valid: the keys should be
        strings, and the values can be integers, strings, floating point values
        or boolean values.
        """
        cdef int intValue = 0
        cdef cbool boolValue = False
        cdef double doubleValue = 0.0
        cdef string strValue
        cdef string keyValue

        if parameterDict:
            for key in parameterDict:
                if type(key) not in (str,unicode):
                    raise ConfigurationParameters("A key is present in the parameter dictionary that is not a string")

                keyValue = B(key)
                val = parameterDict[key]
                if type(val) is int:
                    intValue = val
                    self.m_pParams.setParameter(keyValue, intValue)
                elif type(val) in (str,unicode):
                    strValue = B(val)
                    self.m_pParams.setParameter(keyValue, strValue)
                elif type(val) is float:
                    doubleValue = val
                    self.m_pParams.setParameter(keyValue, doubleValue)
                elif type(val) is bool:
                    boolValue = val
                    self.m_pParams.setParameter(keyValue, boolValue)
                else:
                    raise ConfigurationParametersException("Key '{}' has unsupported value type '{}'".format(key,type(val)))
            
    def __dealloc__(self):
        del self.m_pParams

    def asDict(self):
        """asDict()

        This returns the key/value pairs stored in this instance as a dictionary.
        """
        cdef vector[string] keys
        cdef configurationparameters.TypedParameter param
        cdef cbool boolValue = False
        cdef double doubleValue = 0.0
        cdef string strValue
        cdef string keyValue

        d = { }
        self.m_pParams.getAllKeys(keys)
        for k in keys:
            kStr = S(k)
            if not self.m_pParams.getParameter(k, param):
                raise ConfigurationParametersException("Unable to get parameter for key '{}'".format(kStr))

            if param.isBoolean():
                d[kStr] = param.getBooleanValue()
            elif param.isInteger():
                d[kStr] = param.getIntegerValue()
            elif param.isReal():
                d[kStr] = param.getRealValue()
            elif param.isString():
                d[kStr] = S(param.getStringValue())
            else:
                raise ConfigurationParametersException("Unknown type for key '{}'".format(kStr))

        return d

    @staticmethod
    def fromBytes(paramBytes):
        """fromBytes(paramBytes)

        This method takes in the binary representation of the key/value pairs in
        `paramBytes`, and stores them internally. It is needed to load the information
        from a C++ version of the class.
        """
        cdef configurationparameters.ConfigurationParameters *pParams = NULL
        cdef array[char] buf = chararrayfrombytes(paramBytes)
        cdef serut.MemorySerializer *m = new serut.MemorySerializer(buf.data.as_voidptr, len(paramBytes), NULL, 0)

        try:
            params = ConfigurationParameters({})
            pParams = ConfigurationParameters._getConfigurationParameters(params)
                
            if not pParams.read(deref(m)):
                raise ConfigurationParametersException(S(pParams.getErrorString()))

            return params
        finally:
            del m

class InversionParametersException(Exception):
    """This exception is raised if something goes wrong with the lens inversion
    parameters."""
    pass

cdef shared_ptr[scalesearchparameters.ScaleSearchParameters] _getMassScaleSearchParameters(massScaleSearchType):
    cdef shared_ptr[scalesearchparameters.ScaleSearchParameters] params

    if massScaleSearchType == "regular":
        params.reset(new scalesearchparameters.ScaleSearchParameters(False))
    elif massScaleSearchType == "wide":
        params.reset(new scalesearchparameters.ScaleSearchParameters(True))
    elif massScaleSearchType == "nosearch":
        params.reset(new scalesearchparameters.ScaleSearchParameters())
    else:
        params.reset(new scalesearchparameters.ScaleSearchParameters(
            massScaleSearchType["startFactor"],
            massScaleSearchType["stopFactor"],
            massScaleSearchType["numIterations"],
            massScaleSearchType["firstIterationSteps"],
            massScaleSearchType["nextIterationSteps"]))

    return params

cdef class LensInversionParametersSinglePlaneCPU(object):
    """An internal representation of the parameters for the lens inversion
    procedure, needed to communicate with the C++ based inversion code."""
    
    cdef lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU *m_pParams

    def __cinit__(self):
        self.m_pParams = NULL

    def __dealloc__(self):
        del self.m_pParams

    # each entry in imageList should be a dict with
    # - "images"
    # - "Ds"
    # - "Dds"
    # - (optionally) "params", itself a dict or list
    # - gridInfoOrBasisFunctions:
    #   either a dictionary { "gridSquares": grid squares, "useWeights": flag, "basisFunction": basis type }
    #   or a list of basisfunction info [ { "lens": lens model, "mass": relevant lensing mass, "center": center }, ... ]
    def __init__(self, maxGen, imageList, gridInfoOrBasisFunctions, Dd, zd, massScale,
                 allowNegativeValues = False, baseLens = None, sheetSearch = "nosheet",
                 fitnessObjectParameters = None, massScaleSearchType = "regular"):
        """__init__(maxGen, imageList, gridInfoOrBasisFunctions, Dd, zd, massScale, allowNegativeValues = False, baseLens = None, sheetSearch = "nosheet", fitnessObjectParameters = None,  massScaleSearchType = "regular")
        
        Creates an instance containing the input for the lens inversion method.

        Parameters:
         - `maxGen`: The maximum number of generations in the genetic algorithm. Intended as a fail-safe
           to prevent the algorithm from running indefinitely.
         - `imageList`: this describes the input images, and should be a list of dictionaries where 
           each dictionary has the following entries:

            - ``images``: an :class:`ImagesData <grale.images.ImagesData>` instance, containing
              information about the extended images, point images or null space for example.
            - ``Ds``: angular diameter distance from source to observer for this entry
            - ``Dds``: angular diameter distance between source and lens for this entry
            - ``params`` (optionally): extra parameters that belong to this entry, can be 
              a list or a dictionary.

         - `gridInfoOrBasisFunctions`: for the inversion method you can either provide
           a grid as input, from which basis functions will be derived, or you can specify
           the basis functions directly. 
           
           In the first case, the input should be a dictionary with the following entries:

            - ``gridSquares``: a list of squares in the grid. Each entry is a dictionary
              with ``center`` and ``size`` keys. This size is used directly as the width
              of the used basis function.
            - ``useWeights`` (optional, defaults to ``False``): by default, the initial
              mass of a basis function is the same for a small or large grid cell. If set
              to ``True``, it will be adjusted based on the size of the grid cell, so that
              smaller cells contain less mass.
            - ``basisFunctions`` (optional, defaults to ``plummer``): can be ``plummer``,
              ``square`` or ``gaussian``

           In the second case, it should be a list of dictionaries, where each dictionary
           describes a basis function and should contain the entries:

            - ``lens``: the lens model for this basis function, interpreted as centered on
              the origin.
            - ``mass``: the relevant lensing mass of this basis function. This will be used
              to adjust weights to lie in the relevant total mass range.
            - ``center``: the center (x,y) coordinates at which the lens model should be
              positioned.

         - `Dd`: the angular diameter distance to the lens
         - `zd`: the redshift of the lens, important when using time delay information
         - `massScale`: a mass scale for the optimization to use. The weights of the
           basis functions will be adjusted to lie in a certain range around this
           scale (the width of the range depends on `massScaleSearchType`)
         - `allowNegativeValues`: by default, only positive weights will be assigned to
           the basis functions used in the optimization. If negative weights are allowed
           as well (can be useful when starting from a certain `baseLens`), this can be
           set to ``True``.
         - `baseLens`: if present, the lens inversion method will look for additions to
           this model. The `allowNegativeValues` argument may be useful in this case.
         - `sheetSearch`: can be ``nosheet``, ``genome`` or a lens model, and indicates if a mass sheet
           basis function should be present in the lens inversion method.
         - `fitnessObjectParameters`: parameters that should be passed to the inversion 
           module that will be used.
         - `massScaleSearchType`: default is ``"regular"``, but if set to ``"wide"``, a 
            wider range of masses will be  probed around the specified `massScale` argument. 
            In typical cases, this is not needed, but may be useful when trying to combine 
            strong and weak lensing information. It can also be set to ``"nosearch"`` to
            disable the search for an appropriate scaling of the basis functions. It is also
            possible to set to a dictionary with the following entries, mainly for testing
            purposes:

              - `startFactor`
              - `stopFactor`
              - `numIterations`
              - `firstIterationSteps`
              - `nextIterationSteps`
        """

        cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] imgVector
        cdef vector[grid.GridSquare] gridSquares
        cdef serut.MemorySerializer *mSer = NULL
        cdef array[char] buf
        cdef shared_ptr[imagesdataextended.ImagesDataExtended] pImgDat
        cdef string errorString
        cdef vector2d.Vector2Dd centerVector
        cdef lensinversionparameterssingleplanecpu.BasisFunctionType basisFunctionType
        cdef lensinversionparameterssingleplanecpu.MassSheetSearchType sheetSearchType
        cdef gravitationallens.GravitationalLens *pBaseLens = NULL
        cdef gravitationallens.GravitationalLens *pBasisLensModel = NULL
        cdef gravitationallens.GravitationalLens *pSheetLensModel = NULL
        cdef shared_ptr[gravitationallens.GravitationalLens] sheetLensModel
        cdef configurationparameters.ConfigurationParameters *pFitnessObjectParameters = NULL
        cdef vector[lensinversionbasislensinfo.LensInversionBasisLensInfo] basisLensInfo
        cdef shared_ptr[scalesearchparameters.ScaleSearchParameters] scaleSearchParams

        try:
            # Build the list of extended images
            for imgInfo in imageList:
                imgData = images.ImagesDataExtended(imgInfo["images"])
                imgData.setDs(imgInfo["Ds"])
                imgData.setDds(imgInfo["Dds"])
                
                if "params" in imgInfo:
                    params = imgInfo["params"]
                    paramsSet = False
                    # try a list first, for backward compatibility
                    try:
                        for i in range(len(params)):
                            imgData.setExtraParameter(str(i), params[i])
                        paramsSet = True
                    except:
                        pass

                    if not paramsSet:
                        # Try it as a dictionary
                        for key in params:
                            imgData.setExtraParameter(key, params[key])
        
                # Now that we have the Python ImagesDataExtended object, use
                # toBytes and MemorySerializer to recreate the C++ object
                imgBytes = imgData.toBytes()
                buf = chararrayfrombytes(imgBytes)
                mSer = new serut.MemorySerializer(buf.data.as_voidptr, len(imgBytes), NULL, 0)

                pImgDat.reset(new imagesdataextended.ImagesDataExtended())
                if not pImgDat.get().read(deref(mSer)):
                    errorString = pImgDat.get().getErrorString()
                    raise InversionParametersException(S(errorString))

                imgVector.push_back(pImgDat)    

                del mSer
                mSer = NULL

            # Load a base lens if specified
            if baseLens:
                lensBytes = baseLens.toBytes()
                buf = chararrayfrombytes(lensBytes)
                mSer = new serut.MemorySerializer(buf.data.as_voidptr, len(lensBytes), NULL, 0)

                if not gravitationallens.GravitationalLens.read(deref(mSer), cython.address(pBaseLens), errorString):
                    del mSer
                    raise InversionParametersException(S(errorString))

                del mSer
                mSer = NULL

            # Check the sheet search type
            if sheetSearch == "nosheet":
                sheetSearchType = lensinversionparameterssingleplanecpu.NoSheet
            elif sheetSearch == "genome":
                sheetSearchType = lensinversionparameterssingleplanecpu.Genome
                sheetLensModel = lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU.createDefaultSheetLens(sheetSearchType, Dd)
            elif isinstance(sheetSearch, lenses.GravitationalLens):
                lensBytes = sheetSearch.toBytes()
                buf = chararrayfrombytes(lensBytes)
                mSer = new serut.MemorySerializer(buf.data.as_voidptr, len(lensBytes), NULL, 0)

                if not gravitationallens.GravitationalLens.read(deref(mSer), cython.address(pSheetLensModel), errorString):
                    del mSer
                    raise InversionParametersException(S(errorString))

                del mSer
                mSer = NULL
                sheetLensModel.reset(pSheetLensModel)
            else:
                raise InversionParametersException("Unknown sheet search type '{}', should be 'nosheet' or 'genome', or a lens model".format(sheetSearch))

            if fitnessObjectParameters:
                fitnessObjectParametersObj = ConfigurationParameters(fitnessObjectParameters)
                pFitnessObjectParameters = ConfigurationParameters._getConfigurationParameters(fitnessObjectParametersObj)

            scaleSearchParams = _getMassScaleSearchParameters(massScaleSearchType)

            if type(gridInfoOrBasisFunctions) == dict:
                gridSquareList = gridInfoOrBasisFunctions["gridSquares"]
                basisFunction = gridInfoOrBasisFunctions["basisFunction"] if "basisFunction" in gridInfoOrBasisFunctions else "plummer"
                useWeights = gridInfoOrBasisFunctions["useWeights"] if "useWeights" in gridInfoOrBasisFunctions else False

                # Build the grid square vector
                for square in gridSquareList:
                    cx, cy = square["center"]
                    s = square["size"]

                    centerVector = vector2d.Vector2Dd(float(cx), float(cy))
                    gridSquares.push_back(grid.GridSquare(centerVector, float(s)))

                # Check the basis function type
                if basisFunction == "plummer": 
                    basisFunctionType = lensinversionparameterssingleplanecpu.PlummerBasis
                elif basisFunction == "gaussian": 
                    basisFunctionType = lensinversionparameterssingleplanecpu.GaussBasis
                elif basisFunction == "square":
                    basisFunctionType = lensinversionparameterssingleplanecpu.SquareBasis
                else:
                    raise InversionParametersException("Unknown basis function type '{}', should be 'plummer', 'gaussian' or 'square'".format(basisFunction))

                self.m_pParams = new lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU(maxGen, imgVector, gridSquares,
                                                Dd, zd, massScale, useWeights, basisFunctionType, allowNegativeValues,
                                                pBaseLens, sheetLensModel.get(), pFitnessObjectParameters, deref(scaleSearchParams.get()))

            elif type(gridInfoOrBasisFunctions) == list:

                for entry in gridInfoOrBasisFunctions:
                    cx, cy = entry["center"]
                    relevantLensingMass = entry["mass"]

                    lensBytes = entry["lens"].toBytes()
                    buf = chararrayfrombytes(lensBytes)
                    mSer = new serut.MemorySerializer(buf.data.as_voidptr, len(lensBytes), NULL, 0)

                    if not gravitationallens.GravitationalLens.read(deref(mSer), cython.address(pBasisLensModel), errorString):
                        del mSer
                        raise InversionParametersException(S(errorString))

                    del mSer
                    mSer = NULL

                    LensInversionParametersSinglePlaneCPU._appendHelper(basisLensInfo, pBasisLensModel, cx, cy, relevantLensingMass)

                self.m_pParams = new lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU(maxGen, imgVector, basisLensInfo,
                                                Dd, zd, massScale, allowNegativeValues, pBaseLens, sheetLensModel.get(), 
                                                pFitnessObjectParameters, deref(scaleSearchParams.get()))
            else:
                raise InversionParametersException("Unsupported type for gridInfoOrBasisFunctions parameter, should be dict or list")

        finally:
            # Clean up
            del mSer
            del pBaseLens

    @staticmethod
    cdef _appendHelper(vector[lensinversionbasislensinfo.LensInversionBasisLensInfo] &basisLensInfo,
                       gravitationallens.GravitationalLens *pLensModel, double cx, double cy, double relevantLensingMass):
        cdef shared_ptr[gravitationallens.GravitationalLens] lensModel
        
        lensModel.reset(pLensModel)
        basisLensInfo.push_back(lensinversionbasislensinfo.LensInversionBasisLensInfo(lensModel, vector2d.Vector2Dd(cx, cy), relevantLensingMass))

    cdef _check(self):
        if self.m_pParams == NULL:
            raise InversionParametersException("Internal error: LensInversionParametersSinglePlaneCPU instance has not been set")

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(bytes b)
        
        Creates an instance of this class based on a binary representation. This is
        the inverse function of :func:`toBytes`.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef serut.MemorySerializer *m = new serut.MemorySerializer(buf.data.as_voidptr, len(b), NULL, 0)
        cdef lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU *pParams = NULL
        
        try:
            pParams = new lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU()
            if not pParams.read(deref(m)):
                errStr = S(pParams.getErrorString())
                del pParams
                raise InversionParametersException(errStr)

            r = LensInversionParametersSinglePlaneCPU(0, [], [], 0, 0, 0)
            r.m_pParams = pParams

            return r
        finally:
            del m

    def toBytes(self):
        """toBytes()
        
        Returns a byte array that stores the settings contained in this instance. This
        could be processed again by :func:`fromBytes`.
        """
        cdef serut.VectorSerializer vSer

        self._check()
        if not self.m_pParams.write(vSer):
            raise InversionParametersException(S(self.m_pParams.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

cdef class GAParameters(object):
    """General parameters for the genetic algorithm."""
    
    cdef gaparameters.GAParameters *m_pParams

    def __cinit__(self):
        self.m_pParams = NULL

    def __dealloc__(self):
        del self.m_pParams

    cdef _check(self):
        if self.m_pParams == NULL:
            raise InversionParametersException("Internal error: GAParameters instance has not been set")

    def __init__(self, selectionpressure = 2.5, elitism = True, alwaysincludebest = True, crossoverrate = 0.9):
        """__init__(selectionpressure = 2.5, elitism = True, alwaysincludebest = True, crossoverrate = 0.9)
        
        Initializes an object representing general parameters for the genetic algorithm.
        For more information about the meaning of the arguments, refer to the 
        `documentation <http://research.edm.uhasselt.be/jori/mogal/documentation/classmogal_1_1GeneticAlgorithmParams.html>`_
        of the library that's used for the genetic algorithm.

        """
        self.m_pParams = new gaparameters.GAParameters(selectionpressure, elitism, alwaysincludebest, crossoverrate)

    def getSettings(self):
        """getSettings()

        Returns the settings as a dictionary.
        """
        self._check()
        r = { }
        r["selectionpressure"] = self.m_pParams.getSelectionPressure()
        r["elitism"] = self.m_pParams.getUseElitism()
        r["alwaysincludebest"] = self.m_pParams.getAlwaysIncludeBest()
        r["crossoverrate"] = self.m_pParams.getCrossOverRate()

        return r

    def toBytes(self):
        """toBytes()

        Returns a binary representation of these parameters. Needed for the
        communication with the C++ bases inversion algorithm.
        """
        cdef serut.VectorSerializer vSer

        self._check()
        if not self.m_pParams.write(vSer):
            raise InversionParametersException(S(self.m_pParams.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

cdef class LensInversionParametersMultiPlaneGPU(object):
    """An internal representation of the parameters for the lens multi-plane inversion
    procedure, needed to communicate with the C++ based inversion code."""
    
    cdef lensinversionparametersmultiplanegpu.LensInversionParametersMultiPlaneGPU *m_pParams

    def __cinit__(self):
        self.m_pParams = NULL

    def __dealloc__(self):
        del self.m_pParams

    def __init__(self, cosmology, basisLensesAndRedshifts, imagesAndRedshifts, massEstimate,
                 sheetSearch, fitnessObjectParameters, maxGen, allowNegativeWeights,
                 massScaleSearchType):
        
        cdef cppcosmology.Cosmology cosm
        cdef vector[double] lensRedshifts
        cdef vector[vector[shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo]]] basisLenses
        cdef vector[shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo]] *curPlaneBasisLenses
        cdef shared_ptr[serut.MemorySerializer] mSer
        cdef array[char] buf
        cdef string errorString
        cdef gravitationallens.GravitationalLens *pBasisLensModel = NULL
        cdef shared_ptr[gravitationallens.GravitationalLens] basisLensModel
        cdef shared_ptr[imagesdataextended.ImagesDataExtended] imgDat
        cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] sourceImages
        cdef cbool useSheet
        cdef configurationparameters.ConfigurationParameters *pFitnessObjectParameters = NULL
        cdef shared_ptr[scalesearchparameters.ScaleSearchParameters] scaleSearchParams

        if type(cosmology) != dict:
            cosmology = cosmology.getParameters()

        cosm = cppcosmology.Cosmology(cosmology["h"], cosmology["Omega_m"], cosmology["Omega_r"],
                                      cosmology["Omega_v"], cosmology["w"])
    
        if len(basisLensesAndRedshifts) == 0:
            raise InversionParametersException("No lens planes were specified")

        prevZ = 0
        for planeInfo in basisLensesAndRedshifts:
            z = planeInfo["z"]
            if z <= prevZ:
                raise InversionParametersException("Lens planes should be ordered with stricly increasing redshift")

            prevZ = z

            plane = planeInfo["lenses"]
            if len(plane) == 0:
                raise InversionParametersException("No basis functions in lens plane")

            basisLenses.push_back(vector[shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo]]())
            curPlaneBasisLenses = cython.address(basisLenses.back())

            for entry in plane:
                cx, cy = entry["center"]
                relevantLensingMass = entry["mass"]

                lensBytes = entry["lens"].toBytes()
                buf = chararrayfrombytes(lensBytes)
                mSer.reset(new serut.MemorySerializer(buf.data.as_voidptr, len(lensBytes), NULL, 0))

                if not gravitationallens.GravitationalLens.read(deref(mSer.get()), cython.address(pBasisLensModel), errorString):
                    raise InversionParametersException(S(errorString))
                basisLensModel.reset(pBasisLensModel)

                curPlaneBasisLenses.push_back(
                    shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo](
                        new lensinversionbasislensinfo.LensInversionBasisLensInfo(
                            basisLensModel,
                            vector2d.Vector2Dd(cx, cy), relevantLensingMass)))

        if len(imagesAndRedshifts) == 0:
            raise InversionParametersException("No images specified")
        
        for entry in imagesAndRedshifts:
            img = images.ImagesDataExtended(entry["imgdata"])
            img.setDs(0)  # This distances will be set by the inversion
            img.setDds(0) # This needs to be set to 0
            img.setExtraParameter(entry["z"])

            imgBytes = img.toBytes()
            buf = chararrayfrombytes(imgBytes)
            mSer.reset(new serut.MemorySerializer(buf.data.as_voidptr, len(imgBytes), NULL, 0))

            imgDat.reset(new imagesdataextended.ImagesDataExtended())
            if not imgDat.get().read(deref(mSer)):
                errorString = imgDat.get().getErrorString()
                raise InversionParametersException(S(errorString))

            sourceImages.push_back(imgDat)

        if sheetSearch == "genome":
            useSheet = True
        elif sheetSearch == "nosheet":
            useSheet = False
        else:
            raise InversionParametersException("Bad 'sheetSearch' value, should be 'genome' or 'nosheet'")
        
        if fitnessObjectParameters:
            fitnessObjectParametersObj = ConfigurationParameters(fitnessObjectParameters)
            pFitnessObjectParameters = ConfigurationParameters._getConfigurationParameters(fitnessObjectParametersObj)

        scaleSearchParams = _getMassScaleSearchParameters(massScaleSearchType)

        self.m_pParams = new lensinversionparametersmultiplanegpu.LensInversionParametersMultiPlaneGPU(
            cosm,
            lensRedshifts,
            basisLenses,
            sourceImages,
            massEstimate,
            useSheet,
            pFitnessObjectParameters,
            maxGen,
            allowNegativeWeights,
            deref(scaleSearchParams.get()))