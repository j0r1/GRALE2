"""This module contains classes which will contain various settings
for the lens inversion method. They form the bridge between the Python
side of the code, and the C++ one.
"""
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr, make_shared, unique_ptr, make_unique
from libcpp cimport bool as cbool
from libc.stdint cimport uint64_t
from libcpp.cast cimport dynamic_cast
from libcpp.limits cimport numeric_limits
from libcpp.cmath cimport isnan
from libcpp.utility cimport move
from libcpp.pair cimport pair
from cpython.array cimport array,clone
import cython
from cython.operator cimport dereference as deref
from . import images
from . import lenses
from . import cosmology

import numpy as np
cimport numpy as np

cimport grale.configurationparameters as configurationparameters
cimport grale.lensinversionparameterssingleplanecpu as lensinversionparameterssingleplanecpu
cimport grale.lensinversionparametersmultiplanegpu as lensinversionparametersmultiplanegpu
cimport grale.lensinversionparametersparametricsingleplane as lensinversionparametersparametricsingleplane
cimport grale.serut as serut
cimport grale.imagesdataextended as imagesdataextended
cimport grale.cppgrid as grid
cimport grale.vector2d as vector2d
cimport grale.gravitationallens as gravitationallens
cimport grale.configurationparameters as configurationparameters
cimport grale.eaparameters as eaparameters
cimport grale.lensinversionbasislensinfo as lensinversionbasislensinfo
cimport grale.scalesearchparameters as scalesearchparameters
cimport grale.cppcosmology as cppcosmology
cimport grale.lensgaconverenceparameters as lensgaconverenceparameters
cimport grale.lensgamultipopulationparameters as lensgamultipopulationparameters
cimport grale.errut as errut

include "stringwrappers.pyx"

class ConfigurationParametersException(Exception):
    """This exception is raised in case somethings goes wrong in the
    :class:`ConfigurationParameters` class."""
    pass

cdef class ConfigurationParameters(object):
    """This class describes configuration parameters for a specific lens inversion
    module. It is the Python link to the similarly named C++ class."""

    cdef unique_ptr[configurationparameters.ConfigurationParameters] m_pParams

    def __cinit__(self):
        self.m_pParams = make_unique[configurationparameters.ConfigurationParameters]()

    cdef configurationparameters.ConfigurationParameters * _params(self):
        return self.m_pParams.get()

    @staticmethod
    cdef configurationparameters.ConfigurationParameters * _getConfigurationParameters(ConfigurationParameters confParams):
        return confParams.m_pParams.get()

    def __init__(self, parameterDict = None):
        """__init__(parameterDict = None)

        Initializes this object, storing the key/value entries from `parameterDict`
        (if present). Note that not all dictionaries are valid: the keys should be
        strings, and the values can be integers, strings, floating point values
        or boolean values, or 1D arrays thereof.
        """
        cdef int intValue = 0
        cdef cbool boolValue = False
        cdef double doubleValue = 0.0
        cdef string strValue
        cdef string keyValue
        cdef vector[cbool] boolValues
        cdef vector[int] intValues
        cdef vector[double] doubleValues
        cdef vector[string] stringValues
        cdef np.ndarray[double,ndim=1] realArray
        cdef np.ndarray[np.int32_t,ndim=1] intArray
        cdef np.ndarray[cbool,ndim=1] boolArray
        cdef int idx

        def isArr(x):
            try:
                if type(x) in (str, unicode):
                    return False
                len(x)
                return True
            except:
                return False

        if parameterDict:
            for key in parameterDict:
                if type(key) not in (str,unicode):
                    raise ConfigurationParametersException("A key is present in the parameter dictionary that is not a string")

                keyValue = B(key)
                val = parameterDict[key]
                if val is None:
                    self._params().setParameterEmpty(keyValue)
                    continue

                if not isArr(val):
                    if type(val) is int:
                        intValue = val
                        self._params().setParameter(keyValue, intValue)
                    elif type(val) in (str,unicode):
                        strValue = B(val)
                        self._params().setParameter(keyValue, strValue)
                    elif type(val) is float:
                        doubleValue = val
                        self._params().setParameter(keyValue, doubleValue)
                    elif type(val) is bool:
                        boolValue = val
                        self._params().setParameter(keyValue, boolValue)
                    else:
                        raise ConfigurationParametersException("Key '{}' has unsupported value type '{}'".format(key,type(val)))
                
                else: # Array
                    if type(val) != np.ndarray: # assume it's a list of strings
                        stringValues.clear()
                        for i in val:
                            if type(i) is not str:
                                raise ConfigurationParametersException("Expecting all list elements to be of type string but '{}' is of type '{}'".format(i, type(i)))
    
                            stringValues.push_back(B(i))
                        self._params().setParameter(keyValue, stringValues)
                    else: # Check type of array
                        if len(val.shape) != 1:
                            raise ConfigurationParametersException("Expecting a one dimensional array, but found {} dimensions".format(len(val.shape)))

                        if val.dtype == np.int32:
                            intArray = val
                            intValues.resize(intArray.shape[0])
                            for idx in range(intArray.shape[0]):
                                intValues[idx] = intArray[idx]
                            self._params().setParameter(keyValue, intValues)

                        elif val.dtype == np.float64:
                            realArray = val
                            doubleValues.resize(realArray.shape[0])
                            for idx in range(realArray.shape[0]):
                                doubleValues[idx] = realArray[idx]                                
                            self._params().setParameter(keyValue, doubleValues)

                        elif  val.dtype == np.bool:
                            boolArray = val
                            boolValues.resize(boolArray.shape[0])
                            for idx in range(boolArray.shape[0]):
                                boolValues[idx] = boolArray[idx]
                            self._params().setParameter(keyValue, boolValues)
                        else:                        
                            raise ConfigurationParametersException("Unknown array type, expecting bool, int32 or float64, but got {}".format(val.dtype))
    
    def asDict(self):
        """asDict()

        This returns the key/value pairs stored in this instance as a dictionary.
        """
        cdef vector[string] keys
        cdef const configurationparameters.TypedParameter *param = NULL
        cdef np.ndarray[double,ndim=1] realValues
        cdef np.ndarray[np.int32_t,ndim=1] intValues
        cdef np.ndarray[cbool,ndim=1] boolValues
        cdef const vector[string] *strVec = NULL
        cdef const vector[double] *realVec = NULL
        cdef const vector[int] *intVec = NULL
        cdef const vector[cbool] *boolVec = NULL
        cdef int i

        d = { }
        self._params().getAllKeys(keys)
        for k in keys:
            kStr = S(k)
            param = self._params().getParameter(k)
            if param == NULL:
                raise ConfigurationParametersException("Unable to get parameter for key '{}'".format(kStr))

            # param.dump()
            if not param.isArray():
                if param.isBoolean():
                    d[kStr] = param.getBooleanValue()
                elif param.isInteger():
                    d[kStr] = param.getIntegerValue()
                elif param.isReal():
                    d[kStr] = param.getRealValue()
                elif param.isString():
                    d[kStr] = S(param.getStringValue())
                elif param.isEmpty():
                    d[kStr] = None
                else:
                    raise ConfigurationParametersException("Unknown type for key '{}'".format(kStr))
            else: # Array
                if param.isBoolean():
                    boolVec = cython.address(param.getBooleanValues())
                    boolValues = np.empty([boolVec.size()], dtype=bool)
                    for i in range(boolVec.size()):
                        boolValues[i] = deref(boolVec)[i]
                    d[kStr] = boolValues

                elif param.isInteger():
                    intVec = cython.address(param.getIntegerValues())
                    intValues = np.empty([intVec.size()], dtype=np.int32)
                    for i in range(intVec.size()):
                        intValues[i] = deref(intVec)[i]
                    d[kStr] = intValues

                elif param.isReal():
                    realVec = cython.address(param.getRealValues())
                    realValues = np.empty([realVec.size()], dtype=np.double)
                    for i in range(realVec.size()):
                        realValues[i] = deref(realVec)[i]
                    d[kStr] = realValues

                elif param.isString():
                    values = []
                    strVec = cython.address(param.getStringValues())
                    for i in range(strVec.size()):
                        values.append(S(deref(strVec)[i]))
                    d[kStr] = values

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
        cdef unique_ptr[serut.MemorySerializer] m = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(paramBytes), <void*>NULL, 0)

        params = ConfigurationParameters({})
        pParams = ConfigurationParameters._getConfigurationParameters(params)
            
        if not pParams.read(deref(m)):
            raise ConfigurationParametersException(S(pParams.getErrorString()))

        return params

class InversionParametersException(Exception):
    """This exception is raised if something goes wrong with the lens inversion
    parameters."""
    pass

cdef scalesearchparameters.ScaleSearchParameters* _getMassScaleSearchParameters(massScaleSearchType) except NULL:

    if massScaleSearchType == "regular":
        return new scalesearchparameters.ScaleSearchParameters(False)

    if massScaleSearchType == "wide":
        return new scalesearchparameters.ScaleSearchParameters(True)

    if massScaleSearchType == "nosearch":
        return new scalesearchparameters.ScaleSearchParameters()

    if type(massScaleSearchType) != dict:
        raise InversionParametersException("'massScaleSearchType' must be either 'regular', 'wide', 'nosearch' or a dictionary")

    return new scalesearchparameters.ScaleSearchParameters(
        massScaleSearchType["startFactor"],
        massScaleSearchType["stopFactor"],
        massScaleSearchType["numIterations"],
        massScaleSearchType["firstIterationSteps"],
        massScaleSearchType["nextIterationSteps"])

cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] _createImageVectorFromSinglePlaneImageList(imageList):
    cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] imgVector
    cdef unique_ptr[serut.MemorySerializer] mSer
    cdef array[char] buf
    cdef shared_ptr[imagesdataextended.ImagesDataExtended] pImgDat
    cdef string errorString

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
        mSer = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(imgBytes), <void*>NULL, 0)

        pImgDat.reset(new imagesdataextended.ImagesDataExtended())
        if not pImgDat.get().read(deref(mSer)):
            errorString = pImgDat.get().getErrorString()
            raise InversionParametersException(S(errorString))

        imgVector.push_back(pImgDat)    

    return imgVector

cdef unique_ptr[gravitationallens.GravitationalLens] _createCxxLensFromPyLens(lensModel):
    cdef unique_ptr[gravitationallens.GravitationalLens] cLens
    cdef unique_ptr[serut.MemorySerializer] mSer
    cdef array[char] buf
    cdef string errorString

    if lensModel:
        lensBytes = lensModel.toBytes()
        buf = chararrayfrombytes(lensBytes)
        mSer = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(lensBytes), <void*>NULL, 0)

        if not gravitationallens.GravitationalLens.read(deref(mSer), cLens, errorString):
                raise InversionParametersException(S(errorString))

    return move(cLens)

cdef class LensInversionParametersSinglePlaneCPU(object):
    """An internal representation of the parameters for the lens inversion
    procedure, needed to communicate with the C++ based inversion code."""
    
    cdef unique_ptr[lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU] m_pParams

    # each entry in imageList should be a dict with
    # - "images"
    # - "Ds"
    # - "Dds"
    # - (optionally) "params", itself a dict or list
    # - gridInfoOrBasisFunctions:
    #   either a dictionary { "gridSquares": grid squares, "useWeights": flag, "basisFunction": basis type }
    #   or a list of basisfunction info [ { "lens": lens model, "mass": relevant lensing mass, "center": center }, ... ]
    def __init__(self, imageList, gridInfoOrBasisFunctions, Dd, zd, massScale,
                 allowNegativeValues = False, baseLens = None, sheetSearch = "nosheet",
                 fitnessObjectParameters = None, massScaleSearchType = "regular",
                 randomizeImagePositions = False, initialUncertSeed = 0 
                 ):
        """__init__(maxGen, imageList, gridInfoOrBasisFunctions, Dd, zd, massScale, allowNegativeValues = False, baseLens = None, sheetSearch = "nosheet", fitnessObjectParameters = None,  massScaleSearchType = "regular")
        
        Creates an instance containing the input for the lens inversion method.

        Parameters:
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

        cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] imgVector = _createImageVectorFromSinglePlaneImageList(imageList)
        cdef vector[grid.GridSquare] gridSquares
        cdef unique_ptr[serut.MemorySerializer] mSer
        cdef array[char] buf
        cdef shared_ptr[imagesdataextended.ImagesDataExtended] pImgDat
        cdef string errorString
        cdef vector2d.Vector2Dd centerVector
        cdef lensinversionparameterssingleplanecpu.BasisFunctionType basisFunctionType
        cdef lensinversionparameterssingleplanecpu.MassSheetSearchType sheetSearchType
        cdef unique_ptr[gravitationallens.GravitationalLens] pBaseLens
        cdef unique_ptr[gravitationallens.GravitationalLens] pBasisLensModel
        cdef unique_ptr[gravitationallens.GravitationalLens] pSheetLensModel
        cdef shared_ptr[gravitationallens.GravitationalLens] sheetLensModel
        cdef configurationparameters.ConfigurationParameters *pFitnessObjectParameters = NULL
        cdef vector[lensinversionbasislensinfo.LensInversionBasisLensInfo] basisLensInfo
        cdef shared_ptr[scalesearchparameters.ScaleSearchParameters] scaleSearchParams
        cdef cbool cRandImgPos = randomizeImagePositions
        cdef uint64_t cInitUncertSeed = initialUncertSeed

        # TODO: use _createCxxLensFromPyLens

        # Load a base lens if specified
        if baseLens:
            lensBytes = baseLens.toBytes()
            buf = chararrayfrombytes(lensBytes)
            mSer = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(lensBytes), <void*>NULL, 0)

            if not gravitationallens.GravitationalLens.read(deref(mSer), pBaseLens, errorString):
                raise InversionParametersException(S(errorString))

        # Check the sheet search type
        if sheetSearch == "nosheet":
            sheetSearchType = lensinversionparameterssingleplanecpu.NoSheet
        elif sheetSearch == "genome":
            sheetSearchType = lensinversionparameterssingleplanecpu.Genome
            sheetLensModel = lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU.createDefaultSheetLens(sheetSearchType, Dd)
        elif isinstance(sheetSearch, lenses.GravitationalLens):
            lensBytes = sheetSearch.toBytes()
            buf = chararrayfrombytes(lensBytes)
            mSer = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(lensBytes), <void*>NULL, 0)

            if not gravitationallens.GravitationalLens.read(deref(mSer), pSheetLensModel, errorString):
                raise InversionParametersException(S(errorString))

            sheetLensModel.reset(pSheetLensModel.release())
        else:
            raise InversionParametersException("Unknown sheet search type '{}', should be 'nosheet' or 'genome', or a lens model".format(sheetSearch))

        if fitnessObjectParameters:
            fitnessObjectParametersObj = ConfigurationParameters(fitnessObjectParameters)
            pFitnessObjectParameters = ConfigurationParameters._getConfigurationParameters(fitnessObjectParametersObj)

        scaleSearchParams.reset(_getMassScaleSearchParameters(massScaleSearchType))

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

            self.m_pParams = unique_ptr[lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU](
                    new lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU(imgVector, gridSquares,
                                            Dd, zd, massScale, useWeights, basisFunctionType, allowNegativeValues,
                                            pBaseLens.get(), sheetLensModel.get(), pFitnessObjectParameters, deref(scaleSearchParams.get()),
                                            cRandImgPos, cInitUncertSeed))

        elif type(gridInfoOrBasisFunctions) == list:

            for entry in gridInfoOrBasisFunctions:
                cx, cy = entry["center"]
                relevantLensingMass = entry["mass"]

                lensBytes = entry["lens"].toBytes()
                buf = chararrayfrombytes(lensBytes)
                mSer = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(lensBytes), <void*>NULL, 0)

                if not gravitationallens.GravitationalLens.read(deref(mSer), pBasisLensModel, errorString):
                    raise InversionParametersException(S(errorString))

                LensInversionParametersSinglePlaneCPU._appendHelper(basisLensInfo, pBasisLensModel, cx, cy, relevantLensingMass)

            self.m_pParams = unique_ptr[lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU](
                    new lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU(imgVector, basisLensInfo,
                                            Dd, zd, massScale, allowNegativeValues, pBaseLens.get(), sheetLensModel.get(), 
                                            pFitnessObjectParameters, deref(scaleSearchParams.get()),
                                            cRandImgPos, cInitUncertSeed))
        else:
            raise InversionParametersException("Unsupported type for gridInfoOrBasisFunctions parameter, should be dict or list")

    @staticmethod
    cdef _appendHelper(vector[lensinversionbasislensinfo.LensInversionBasisLensInfo] &basisLensInfo,
                       unique_ptr[gravitationallens.GravitationalLens] &pLensModel, double cx, double cy, double relevantLensingMass):
        cdef shared_ptr[gravitationallens.GravitationalLens] lensModel
        
        lensModel.reset(pLensModel.release())
        basisLensInfo.push_back(lensinversionbasislensinfo.LensInversionBasisLensInfo(lensModel, vector2d.Vector2Dd(cx, cy), relevantLensingMass))

    cdef _check(self):
        if self.m_pParams.get() == NULL:
            raise InversionParametersException("Internal error: LensInversionParametersSinglePlaneCPU instance has not been set")

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(bytes b)
        
        Creates an instance of this class based on a binary representation. This is
        the inverse function of :func:`toBytes`.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef unique_ptr[serut.MemorySerializer] m = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(b), <void*>NULL, 0)
        cdef unique_ptr[lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU] pParams = make_unique[lensinversionparameterssingleplanecpu.LensInversionParametersSinglePlaneCPU]()
        
        if not deref(pParams).read(deref(m)):
            raise InversionParametersException(S(deref(pParams).getErrorString()))

        r = LensInversionParametersSinglePlaneCPU(0, [], [], 0, 0, 0)
        r.m_pParams.swap(pParams)

        return r

    def toBytes(self):
        """toBytes()
        
        Returns a byte array that stores the settings contained in this instance. This
        could be processed again by :func:`fromBytes`.
        """
        cdef serut.VectorSerializer vSer

        self._check()
        if not deref(self.m_pParams).write(vSer):
            raise InversionParametersException(S(deref(self.m_pParams).getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

# TODO: move some code to common base class

cdef class EAParameters(object):
    """Base class for EA algorithm parameters."""

    cdef unique_ptr[eaparameters.EAParameters] m_pParams

    cdef _check(self):
        if self.m_pParams.get() == NULL:
            raise InversionParametersException("Internal error: parameters instance has not been set")

    def getSettings(self):
        """getSettings()

        Returns the settings as a dictionary.
        """
        self._check()
        r = { }
        self._fillInSettings(r)

        return r

    def _fillInSettings(self, r):
        raise Exception("Internal error: should be overridden in derived class")

    def toBytes(self):
        """toBytes()

        Returns a binary representation of these parameters. Needed for the
        communication with the C++ bases inversion algorithm.
        """
        cdef serut.VectorSerializer vSer
        cdef errut.bool_t r;

        self._check()
        r = deref(self.m_pParams).write(vSer);
        if not r.success():
            raise InversionParametersException(S(r.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

cdef class RNDParameters(EAParameters):
    """TODO"""
    
    def __init__(self, scale = 1e-5):
        """__init__()
        
        TODO
        """
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.RNDParameters(scale))

    def _fillInSettings(self, r):
        cdef eaparameters.RNDParametersPtrConst pParams = dynamic_cast[eaparameters.RNDParametersPtrConst](self.m_pParams.get())
        if not pParams:
            raise Exception("Internal error: can't dynamic_cast parameters to correct type")

        r["scale"] = deref(pParams).getScale()

cdef class EATestParameters(EAParameters):
    """General parameters for the TEST algorithm."""

    def __init__(self):
        """__init__()

        TODO
        """
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.EATestParameters())

    def _fillInSettings(self, r):
        cdef eaparameters.EATestParametersPtrConst pParams = dynamic_cast[eaparameters.EATestParametersPtrConst](self.m_pParams.get())
        if not pParams:
            raise Exception("Internal error: can't dynamic_cast parameters to correct type")

cdef class DEParameters(EAParameters):
    """General parameters for the DE algorithm."""
    
    def __init__(self, F = 0.5, CR = 0.5, needstrictlybetter = True): # TODO: what are good defaults?
        """__init__()
        
        TODO
        """
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.DEParameters(F, CR, needstrictlybetter))

    def _fillInSettings(self, r):
        cdef eaparameters.DEParametersPtrConst pParams = dynamic_cast[eaparameters.DEParametersPtrConst](self.m_pParams.get())
        if not pParams:
            raise Exception("Internal error: can't dynamic_cast parameters to correct type")

        r["F"] = deref(pParams).getF()
        r["CR"] = deref(pParams).getCR()
        r["needstrictlybetter"] = deref(pParams).getNeedStrictlyBetter()

cdef class JADEParameters(EAParameters):
    """General parameters for the JADE algorithm."""
    
    def __init__(self, p = 0.05, c = 0.1, usearchive = True, initmeanF = 0.5, initmeanCR = 0.5, needstrictlybetter = True):
        """__init__()
        
        TODO
        """
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.JADEParameters(p, c, usearchive, initmeanF, initmeanCR, needstrictlybetter))

    def _fillInSettings(self, r):
        cdef eaparameters.JADEParametersPtrConst pParams = dynamic_cast[eaparameters.JADEParametersPtrConst](self.m_pParams.get())
        if not pParams:
            raise Exception("Internal error: can't dynamic_cast parameters to correct type")

        r["p"] = deref(pParams).getBestFraction_p()
        r["c"] = deref(pParams).getParameterUpdateFraction_c()
        r["usearchive"] = deref(pParams).useExternalArchive()
        r["initmeanF"] = deref(pParams).getInitialMeanF()
        r["initmeanCR"] = deref(pParams).getInitialMeanCR()
        r["needstrictlybetter"] = deref(pParams).getNeedStrictlyBetter()

cdef class GAParameters(EAParameters):
    """General parameters for the genetic algorithm."""

    def __init__(self, selectionpressure = 2.5, elitism = True, alwaysincludebest = True, crossoverrate = 0.9, smallmutationsize = -1):
        """__init__(selectionpressure = 2.5, elitism = True, alwaysincludebest = True, crossoverrate = 0.9, smallmutationsize = -1)
        
        Initializes an object representing general parameters for the genetic algorithm.
        For more information about the meaning of the arguments, refer to the 
        `documentation <http://research.edm.uhasselt.be/jori/mogal/documentation/classmogal_1_1GeneticAlgorithmParams.html>`_
        of the library that's used for the genetic algorithm.

        TODO: update this documentation
        """
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.GAParameters(selectionpressure, elitism, alwaysincludebest, crossoverrate, smallmutationsize))

    def _fillInSettings(self, r):
        cdef eaparameters.GAParametersPtrConst pParams = dynamic_cast[eaparameters.GAParametersPtrConst](self.m_pParams.get())
        if not pParams:
            raise Exception("Internal error: can't dynamic_cast parameters to correct type")

        r["selectionpressure"] = deref(pParams).getSelectionPressure()
        r["elitism"] = deref(pParams).getUseElitism()
        r["alwaysincludebest"] = deref(pParams).getAlwaysIncludeBest()
        r["crossoverrate"] = deref(pParams).getCrossOverRate()
        r["smallmutationsize"] = deref(pParams).getSmallMutationSize()

cdef class NSGA2Parameters(EAParameters):
    """TODO"""

    def __init__(self, smallmutationsize = -1):
        """__init__(smallmutationsize = -1)

        TODO
        """
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.NSGA2Parameters(smallmutationsize))

    def _fillInSettings(self, r):
        cdef eaparameters.NSGA2ParametersPtrConst pParams = dynamic_cast[eaparameters.NSGA2ParametersPtrConst](self.m_pParams.get())
        if not pParams:
            raise Exception("Internal error: can't dynamic_cast parameters to correct type")

        r["smallmutationsize"] = deref(pParams).getSmallMutationSize()

cdef class NSGA2DELikeCrossoverParameters(EAParameters):
    """TODO"""

    def __init__(self, useextraparent = True, F = "random", CR = "random"):
        """__init__()

        TODO
        """
        cdef float cF, cCR;

        cF = numeric_limits[float].quiet_NaN() if F == "random" else F
        cCR = numeric_limits[float].quiet_NaN() if CR == "random" else CR
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.NSGA2DELikeCrossoverParameters(useextraparent, cF, cCR))

    def _fillInSettings(self, r):
        cdef eaparameters.NSGA2DELikeCrossoverParametersPtrConst pParams = dynamic_cast[eaparameters.NSGA2DELikeCrossoverParametersPtrConst](self.m_pParams.get())
        cdef float F, CR
        if not pParams:
            raise Exception("Internal error: can't dynamic_cast parameters to correct type")

        F = deref(pParams).getF()
        CR = deref(pParams).getCR()

        r["useextraparent"] = deref(pParams).useExtraParent()
        r["F"] = "random" if isnan(F) else F
        r["CR"] = "random" if isnan(CR) else CR

cdef class MCMCParameters(EAParameters):
    """TODO"""
    def __init__(self, a = 2.0, samplesfilename = "", samplegenerations = 0, burningenerations = 0,
                 annealgenerationsscale = 0, alpha0 = 0.1, alphamax = 1.0):
        """__init()__"

        TODO
        """
        cdef double cA, cAlpha0, cAlphaMax
        cdef string cFileName
        cdef size_t cGen, cAnnealGen, cBurnGen
        
        cA = a
        cAlpha0 = alpha0
        cAlphaMax = alphamax
        cFileName = B(samplesfilename)
        cGen = samplegenerations
        cBurnGen = burningenerations
        cAnnealGen = annealgenerationsscale
        self.m_pParams = unique_ptr[eaparameters.EAParameters](new eaparameters.MCMCParameters(cA, cFileName, cGen, cBurnGen, cAnnealGen, cAlpha0, cAlphaMax))

    def _fillInSettings(self, r):
        cdef eaparameters.MCMCParametersPtrConst pParams = dynamic_cast[eaparameters.MCMCParametersPtrConst](self.m_pParams.get())
        r["a"] = deref(pParams).getGoodmanWeare_a()
        r["samplesfilename"] = S(deref(pParams).getSamplesFilename())
        r["samplegenerations"] = deref(pParams).getSampleGenerations()
        r["burningenerations"] = deref(pParams).getBurnInGenerations()
        r["annealgenerationsscale"] = deref(pParams).getAnnealGenerationsTimeScale()
        r["alpha0"] = deref(pParams).getAnnealAlpha0()
        r["alphamax"] = deref(pParams).getAnnealAlphaMax()

cdef class LensInversionParametersMultiPlaneGPU(object):
    """An internal representation of the parameters for the lens multi-plane inversion
    procedure, needed to communicate with the C++ based inversion code."""
    
    cdef unique_ptr[lensinversionparametersmultiplanegpu.LensInversionParametersMultiPlaneGPU] m_pParams

    def __init__(self, cosmology = cosmology.Cosmology(0.7, 0.3, 0, 0.7),
                 basisLensesAndRedshifts = [], imagesAndRedshifts = [], baseLensForPlane = [],
                 massEstimate = 0, sheetSearch = "nosheet", fitnessObjectParameters = None,
                 allowNegativeWeights = False, massScaleSearchType = "regular",
                 deviceIndex = "rotate"):
        """__init__(cosmology = cosmology.Cosmology(0.7, 0.3, 0, 0.7), basisLensesAndRedshifts = [], imagesAndRedshifts = [], massEstimate = 0, sheetSearch = "nosheet", fitnessObjectParameters = None, maxGen = 0, allowNegativeWeights = False, massScaleSearchType = "regular", deviceIndex = "rotate")

        Parameters:

         - `cosmology`: the cosmological model to use to transform redshifts into angular diameter
           distances.
         - `basisLensesAndRedshifts`: a list of dictionaries with entries `z` and `lenses`, 
           representing the lens planes and the basis functions therein. Each entry in the `lenses`
           list should be a dictionary with `lens`, `center` and `mass` entries, as with the
           single plane inversion.
         - `imagesAndRedshifts`: a list of dictionaries, each having an entry `z` describing the
           redshift of the images, `images` with an :class:`ImagesData <grale.images.ImagesData>`
           instance, and `params` listing additional parameters (such as a `type`).
         - `baseLensForPlane`: if present, this list should contain one base lens (possibly ``None``)
           for each lens plane.
         - `massEstimate`: a mass scale for the optimization to use. The weights of the
           basis functions will be adjusted to lie in a certain range around this
           scale (the width of the range depends on `massScaleSearchType`)
         - `sheetSearch`: can be ``nosheet`` or ``genome``.
         - `fitnessObjectParameters`: parameters that should be passed to the inversion 
           module that will be used.
         - `allowNegativeWeights`: by default, only positive weights will be assigned to
           the basis functions used in the optimization. If negative weights are allowed
           as well, this can be set to ``True``.
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

         - `deviceIndex`: this multi-plane inversion uses a GPU to back-project the image
           data, and by setting a specific number, a specific device can be specified. To
           allow multiple GPUs to be used automatically, you can leave this to ``"auto"``
           and use an :mod:`inverter <grale.inverters>` with as many processes as you have
           GPUs.
        """
        
        cdef cppcosmology.Cosmology cosm
        cdef vector[double] lensRedshifts
        cdef vector[vector[shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo]]] basisLenses
        cdef vector[shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo]] *curPlaneBasisLenses
        cdef vector[shared_ptr[gravitationallens.GravitationalLens]] baseLensesPerPlane
        cdef shared_ptr[serut.MemorySerializer] mSer
        cdef array[char] buf
        cdef string errorString
        cdef unique_ptr[gravitationallens.GravitationalLens] pBasisLensModel
        cdef shared_ptr[gravitationallens.GravitationalLens] basisLensModel
        cdef shared_ptr[imagesdataextended.ImagesDataExtended] imgDat
        cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] sourceImages
        cdef cbool useSheet
        cdef configurationparameters.ConfigurationParameters *pFitnessObjectParameters = NULL
        cdef shared_ptr[scalesearchparameters.ScaleSearchParameters] scaleSearchParams
        cdef int devIdx

        if type(cosmology) != dict:
            cosmology = cosmology.getParameters()

        cosm = cppcosmology.Cosmology(cosmology["h"], cosmology["Omega_m"], cosmology["Omega_r"],
                                      cosmology["Omega_v"], cosmology["w"])
    
        prevZ = 0
        for planeInfo in basisLensesAndRedshifts:
            z = planeInfo["z"]
            if z <= prevZ:
                raise InversionParametersException("Lens planes should be ordered with stricly increasing redshift")

            prevZ = z
            lensRedshifts.push_back(z)

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

                if not gravitationallens.GravitationalLens.read(deref(mSer.get()), pBasisLensModel, errorString):
                    raise InversionParametersException(S(errorString))
                basisLensModel.reset(pBasisLensModel.release())

                curPlaneBasisLenses.push_back(shared_ptr[lensinversionbasislensinfo.LensInversionBasisLensInfo](
                    new lensinversionbasislensinfo.LensInversionBasisLensInfo(
                            basisLensModel,
                            vector2d.Vector2Dd(cx, cy), relevantLensingMass)))

        if baseLensForPlane:
            if len(baseLensForPlane) != lensRedshifts.size():
                raise InversionParametersException("The number of base lenses should equal the number of lens planes")

            baseLensesPerPlane.resize(lensRedshifts.size()) # Initializes each lens to nullptr
            for i in range(len(baseLensForPlane)):
                l = baseLensForPlane[i]

                if l: # Can be None
                    lensBytes = l.toBytes()
                    buf = chararrayfrombytes(lensBytes)
                    mSer.reset(new serut.MemorySerializer(buf.data.as_voidptr, len(lensBytes), NULL, 0))

                    if not gravitationallens.GravitationalLens.read(deref(mSer.get()), pBasisLensModel, errorString):
                        raise InversionParametersException(S(errorString))

                    baseLensesPerPlane[i].reset(pBasisLensModel.release())

        for entry in imagesAndRedshifts:
            img = images.ImagesDataExtended(entry["images"])
            img.setDs(0)  # This distances will be set by the inversion
            img.setDds(0) # This needs to be set to 0
            if "params" in entry:
                params = entry["params"]
                for key in params:
                    img.setExtraParameter(key, params[key])

            img.setExtraParameter("z", entry["z"])

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

        scaleSearchParams.reset(_getMassScaleSearchParameters(massScaleSearchType))

        if deviceIndex == "rotate":
            devIdx = -1
        else:
            if deviceIndex < 0:
                raise InversionParametersException("Device index can't be negative")
            devIdx = deviceIndex

        self.m_pParams = unique_ptr[lensinversionparametersmultiplanegpu.LensInversionParametersMultiPlaneGPU](
            new lensinversionparametersmultiplanegpu.LensInversionParametersMultiPlaneGPU(
                cosm,
                lensRedshifts,
                basisLenses,
                baseLensesPerPlane,
                sourceImages,
                massEstimate,
                useSheet,
                pFitnessObjectParameters,
                allowNegativeWeights,
                deref(scaleSearchParams.get()),
                devIdx))

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(bytes b)

        Creates an instance of this class based on a binary representation. This is
        the inverse function of :func:`toBytes`.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef unique_ptr[serut.MemorySerializer] m = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(b), <void*>NULL, 0)

        r = LensInversionParametersMultiPlaneGPU()
        if not deref(r.m_pParams).read(deref(m)):
            raise InversionParametersException(S(deref(r.m_pParams).getErrorString()))

        return r

    def toBytes(self):
        """toBytes()

        Returns a byte array that stores the settings contained in this instance. This
        could be processed again by :func:`fromBytes`.
        """
        cdef serut.VectorSerializer vSer

        if not deref(self.m_pParams).write(vSer):
            raise InversionParametersException(S(deref(self.m_pParams).getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

cdef vector[int] _createIntVectorFromList(l):
    cdef vector[int] result
    cdef int cValue
    for x in l:
        cValue = x
        result.push_back(cValue)
    return result

cdef vector[float] _createFloatVectorFromList(l):
    cdef vector[float] result
    cdef float cValue
    for x in l:
        cValue = x
        result.push_back(cValue)
    return result

cdef class LensInversionParametersParametricSinglePlane(object):

    cdef unique_ptr[lensinversionparametersparametricsingleplane.LensInversionParametersParametricSinglePlane] m_pParams

    def __init__(self, inputImages, Dd, zd,
                  templateLens, deflScale, potScale, offsets, initMin, initMax, hardMin, hardMax,
                  infOnBoundsViolation,
                  fitnessObjectParameters, deviceIndex = "rotate",
                  randomizeImagePositions = False, initialUncertSeed = 0,
                  originParametersMap = [], numOriginParams = 0):

        cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] imgVector = _createImageVectorFromSinglePlaneImageList(inputImages)
        cdef double cDd = Dd
        cdef double cZd = zd
        cdef unique_ptr[gravitationallens.GravitationalLens] cTemplateLens = _createCxxLensFromPyLens(templateLens)
        cdef double cDeflScale = deflScale
        cdef double cPotScale = potScale        
        cdef vector[int] cOffsets = _createIntVectorFromList(offsets)
        cdef vector[float] cInitMin = _createFloatVectorFromList(initMin)
        cdef vector[float] cInitMax = _createFloatVectorFromList(initMax)
        cdef vector[float] cHardMin = _createFloatVectorFromList(hardMin)
        cdef vector[float] cHardMax = _createFloatVectorFromList(hardMax)
        cdef configurationparameters.ConfigurationParameters *pFitnessObjectParameters = NULL
        cdef int devIdx
        cdef cbool cInfOnBoundsViolation = infOnBoundsViolation
        cdef cbool cRandomizeInputPos = randomizeImagePositions
        cdef uint64_t cInitialUncertSeed = initialUncertSeed
        cdef vector[pair[size_t,string]] cOriginParamMapping
        cdef size_t cNumOriginParams = numOriginParams
        cdef size_t tmpSize
        cdef string tmpStr

        if inputImages is None:
            return

        fitnessObjectParametersObj = ConfigurationParameters(fitnessObjectParameters)
        pFitnessObjectParameters = ConfigurationParameters._getConfigurationParameters(fitnessObjectParametersObj)

        if deviceIndex == "rotate":
            devIdx = -1
        else:
            if deviceIndex < 0:
                raise InversionParametersException("Device index can't be negative")
            devIdx = deviceIndex

        for entry in originParametersMap:
            pos, code = entry
            if pos < 0:
                raise InversionParametersException("Invalid negative position in origin parameters map")
            tmpSize = pos
            tmpStr = B(code)
            cOriginParamMapping.push_back(pair[size_t,string](tmpSize, tmpStr))

        self.m_pParams = unique_ptr[lensinversionparametersparametricsingleplane.LensInversionParametersParametricSinglePlane](
            new lensinversionparametersparametricsingleplane.LensInversionParametersParametricSinglePlane(
                imgVector, cDd, cZd, deref(cTemplateLens), cDeflScale, cPotScale,
                cOffsets, cInitMin, cInitMax, cHardMin, cHardMax, cInfOnBoundsViolation, deref(pFitnessObjectParameters),
                devIdx, cRandomizeInputPos, cInitialUncertSeed, cOriginParamMapping, cNumOriginParams
            )
        )

    cdef _check(self):
        if self.m_pParams.get() == NULL:
            raise InversionParametersException("Internal error: LensInversionParametersParametricSinglePlane instance has not been set")

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(bytes b)
        
        Creates an instance of this class based on a binary representation. This is
        the inverse function of :func:`toBytes`.
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef unique_ptr[serut.MemorySerializer] m = make_unique[serut.MemorySerializer](buf.data.as_voidptr, len(b), <void*>NULL, 0)
        cdef unique_ptr[lensinversionparametersparametricsingleplane.LensInversionParametersParametricSinglePlane] pParams = make_unique[lensinversionparametersparametricsingleplane.LensInversionParametersParametricSinglePlane]()
        
        if not deref(pParams).read(deref(m)):
            raise InversionParametersException(S(deref(pParams).getErrorString()))

        r = LensInversionParametersParametricSinglePlane(None, None, None, None, None, None, None,
                                                         None, None, None, None, None, None, None)
        r.m_pParams.swap(pParams)

        return r

    def toBytes(self):
        """toBytes()
        
        Returns a byte array that stores the settings contained in this instance. This
        could be processed again by :func:`fromBytes`.
        """
        cdef serut.VectorSerializer vSer

        self._check()
        if not deref(self.m_pParams).write(vSer):
            raise InversionParametersException(S(deref(self.m_pParams).getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

class ConvergenceParametersException(Exception):
    """This exception is raised in case somethings goes wrong in the
    :class:`ConvergenceParameters` class."""
    pass

cdef class ConvergenceParameters(object):
    """This class is used to control the settings that determine when a
    part of the GA/EA will be considered as converged. In the classic
    GA, two such stages are used, one with larger mutations, and one with
    smaller, each also with a different convergence factor.

    The settings are controlled with a dictionary with these entries:

     - ``maximumgenerations``: if the number of generations of the GA exceeds
       this number, it will stop (defaults to 16384).
     - ``historysize``: to determine how the fitness measures are converging,
       the current fitness values are compared to those of a number of generations
       ago. This number is what's described by the ``historysize`` value (defaults
       to 250). If the relative change of all fitness measures (in a multi-objective
       setting) are below a certain threshold, convergence is reached.
     - `convergencefactor`: this contains the threshold mentioned above.

    """

    cdef unique_ptr[lensgaconverenceparameters.LensGAConvergenceParameters] m_params
    cdef lensgaconverenceparameters.LensGAConvergenceParameters *m_pParams

    def __cinit__(self):
        self.m_params = make_unique[lensgaconverenceparameters.LensGAConvergenceParameters]()
        self.m_pParams = self.m_params.get()

    def __init__(self, parameterDict = None, eaType = None): # eaType is only used for the defaults
        if not parameterDict:
            parameterDict = {
                "maximumgenerations": 16384,
                "historysize": 250,
                "convergencefactor": 0.05,
            }

        self.fromDict(parameterDict)

    def fromDict(self, d):
        
        cdef vector[double] factors, sizes;

        knownKeys = [ "maximumgenerations",
                      "historysize",
                      "convergencefactor" ]

        for k in d:
            if not k in knownKeys:
                raise ConvergenceParametersException("Unknown key {}, known are {}".format(k, knownKeys))

        self.m_pParams.setMaximumNumberOfGenerations(d["maximumgenerations"])
        self.m_pParams.setHistorySize(d["historysize"])
        self.m_pParams.setConvergenceFactor(d["convergencefactor"])

    def toDict(self):
        d = { }
        d["maximumgenerations"] = self.m_pParams.getMaximumNumberOfGenerations()
        d["historysize"] = self.m_pParams.getConvergenceHistorySize()
        d["convergencefactor"] = self.m_pParams.getConvergenceFactor()

        return d

    def toBytes(self):
        """toBytes()

        Returns a byte array that stores the settings contained in this instance.
        """
        cdef serut.VectorSerializer vSer

        if not self.m_pParams.write(vSer):
            raise ConvergenceParametersException(S(self.m_pParams.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

class MultiPopulationParametersException(Exception):
    """This exception is raised in case somethings goes wrong in the
    :class:`MultiPopulationParameters` class."""
    pass

cdef class MultiPopulationParameters(object):
    """TODO"""

    cdef unique_ptr[lensgamultipopulationparameters.LensGAMultiPopulationParameters] m_params
    cdef lensgamultipopulationparameters.LensGAMultiPopulationParameters *m_pParams

    def __cinit__(self):
        self.m_params = make_unique[lensgamultipopulationparameters.LensGAMultiPopulationParameters]()
        self.m_pParams = self.m_params.get()

    def __init__(self, parameterDict):
        self.fromDict(parameterDict)

    def fromDict(self, d):
        knownKeys = [ "populations", "skipgenerations", "generationfraction", "migrants" ]
        for k in d:
            if not k in knownKeys:
                raise MultiPopulationParametersException("Key '{}' is not supported".format(k))

        if set(knownKeys) != d.keys():
            raise MultiPopulationParametersException("Not all keys are present, need: '{}'".format(knownKeys))

        self.m_pParams.setNumberOfPopulations(d["populations"])
        self.m_pParams.setNumberOfInitialGenerationsToSkip(d["skipgenerations"])
        self.m_pParams.setMigrationGenerationFraction(d["generationfraction"])
        self.m_pParams.setNumberOfIndividualsToLeavePopulation(d["migrants"])

    def toDict(self):
        d = { }
        d["populations"] = self.m_pParams.getNumberOfPopulations()
        d["skipgenerations"] = self.m_pParams.getNumberOfInitialGenerationsToSkip()
        d["generationfraction"] = self.getMigrationGenerationFraction()
        d["migrants"] = self.getNumberOfIndividualsToLeavePopulation()
        return d

    def toBytes(self):
        """toBytes()

        Returns a byte array that stores the settings contained in this instance.
        """
        cdef serut.VectorSerializer vSer

        if not self.m_pParams.write(vSer):
            raise MultiPopulationParametersException(S(self.m_pParams.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

