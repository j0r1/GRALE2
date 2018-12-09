""" TODO

"""
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr
from libcpp cimport bool as cbool
from cpython.array cimport array,clone
import cython
from cython.operator cimport dereference as deref
from . import images

cimport grale.configurationparameters as configurationparameters
cimport grale.gridlensinversionparameters as gridlensinversionparameters
cimport grale.serut as serut
cimport grale.imagesdataextended as imagesdataextended
cimport grale.cppgrid as grid
cimport grale.vector2d as vector2d
cimport grale.gravitationallens as gravitationallens
cimport grale.configurationparameters as configurationparameters
cimport grale.gaparameters as gaparameters

include "stringwrappers.pyx"

class ConfigurationParametersException(Exception):
    """TODO:"""
    pass

cdef class ConfigurationParameters(object):
    """TODO:"""

    cdef configurationparameters.ConfigurationParameters *m_pParams

    def __cinit__(self):
        self.m_pParams = new configurationparameters.ConfigurationParameters()

    @staticmethod
    cdef configurationparameters.ConfigurationParameters * _getConfigurationParameters(confParams):
        return confParams.m_pParams

    def __init__(self, parameterDict = None):
        """TODO:"""
        cdef int intValue = 0
        cdef cbool boolValue = False
        cdef double doubleValue = 0.0
        cdef string strValue = ""
        cdef string keyValue = ""

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
        """TODO:"""
        cdef vector[string] keys
        cdef configurationparameters.TypedParameter param
        cdef cbool boolValue = False
        cdef double doubleValue = 0.0
        cdef string strValue = ""
        cdef string keyValue = ""

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
        """TODO:"""
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
    """TODO:"""
    pass

cdef class GridLensInversionParameters(object):
    """TODO:"""
    
    cdef gridlensinversionparameters.GridLensInversionParameters *m_pParams

    def __cinit__(self):
        self.m_pParams = NULL

    def __dealloc__(self):
        del self.m_pParams

    # each entry in imageList should be a dict with
    # - "images"
    # - "Ds"
    # - "Dds"
    # - (optionally) "params", itself a dict or list
    def __init__(self, maxGen, imageList, gridSquareList, Dd, zd, massScale, useWeights = False, basisFunction = "plummer",
                 allowNegativeValues = False, baseLens = None, sheetSearch = "nosheet",
                 fitnessObjectParameters = None, wideSearch = False):
        """TODO:"""

        cdef vector[shared_ptr[imagesdataextended.ImagesDataExtended]] imgVector
        cdef vector[grid.GridSquare] gridSquares
        cdef serut.MemorySerializer *mSer
        cdef array[char] buf
        cdef shared_ptr[imagesdataextended.ImagesDataExtended] pImgDat
        cdef string errorString
        cdef vector2d.Vector2Dd centerVector
        cdef gridlensinversionparameters.BasisFunctionType basisFunctionType
        cdef gridlensinversionparameters.MassSheetSearchType sheetSearchType
        cdef gravitationallens.GravitationalLens *pBaseLens = NULL
        cdef configurationparameters.ConfigurationParameters *pFitnessObjectParameters = NULL

        mSer = NULL

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

            # Build the grid square vector
            for square in gridSquareList:
                cx, cy = square["center"]
                s = square["size"]

                centerVector = vector2d.Vector2Dd(float(cx), float(cy))
                gridSquares.push_back(grid.GridSquare(centerVector, float(s)))

            # Check the basis function type
            if basisFunction == "plummer": 
                basisFunctionType = gridlensinversionparameters.PlummerBasis
            elif basisFunction == "gaussian": 
                basisFunctionType = gridlensinversionparameters.GaussBasis
            elif basisFunction == "square":
                basisFunctionType = gridlensinversionparameters.SquareBasis
            else:
                raise InversionParametersException("Unknown basis function type '{}', should be 'plummer', 'gaussian' or 'square'".format(basisFunction))

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
                sheetSearchType = gridlensinversionparameters.NoSheet
            elif sheetSearch == "genome":
                sheetSearchType = gridlensinversionparameters.Genome
            else:
                raise InversionParametersException("Unknown sheet search type '{}', should be 'nosheet' or 'genome'".format(sheetSearch))

            if fitnessObjectParameters:
                fitnessObjectParametersObj = ConfigurationParameters(fitnessObjectParameters)
                pFitnessObjectParameters = ConfigurationParameters._getConfigurationParameters(fitnessObjectParametersObj)

            self.m_pParams = new gridlensinversionparameters.GridLensInversionParameters(maxGen, imgVector, gridSquares,
                                                Dd, zd, massScale, useWeights, basisFunctionType, allowNegativeValues,
                                                pBaseLens, sheetSearchType, pFitnessObjectParameters, wideSearch)


        finally:
            # Clean up
            del mSer
            del pBaseLens

    cdef _check(self):
        if self.m_pParams == NULL:
            raise InversionParametersException("Internal error: GridLensInversionParameters instance has not been set")

    @staticmethod
    def fromBytes(bytes b):
        """fromBytes(bytes b)
        
        TODO
        """
        cdef array[char] buf = chararrayfrombytes(b)
        cdef serut.MemorySerializer *m = new serut.MemorySerializer(buf.data.as_voidptr, len(b), NULL, 0)
        cdef gridlensinversionparameters.GridLensInversionParameters *pParams = NULL
        
        try:
            pParams = new gridlensinversionparameters.GridLensInversionParameters()
            if not pParams.read(deref(m)):
                errStr = S(pParams.getErrorString())
                del pParams
                raise InversionParametersException(errStr)

            r = GridLensInversionParameters()
            r.m_pParams = pParams

            return r
        finally:
            del m

    def toBytes(self):
        """toBytes()
        
        TODO
        """
        cdef serut.VectorSerializer vSer

        self._check()
        if not self.m_pParams.write(vSer):
            raise InversionParametersException(S(self.m_pParams.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

cdef class GAParameters(object):
    """TODO:"""
    
    cdef gaparameters.GAParameters *m_pParams

    def __cinit__(self):
        self.m_pParams = NULL

    def __dealloc__(self):
        del self.m_pParams

    cdef _check(self):
        if self.m_pParams == NULL:
            raise InversionParametersException("Internal error: GAParameters instance has not been set")

    def __init__(self, selectionpressure = 2.5, elitism = True, alwaysincludebest = True, crossoverrate = 0.9):
        """TODO:"""
        self.m_pParams = new gaparameters.GAParameters(selectionpressure, elitism, alwaysincludebest, crossoverrate)

    def getSettings(self):
        """TODO:"""
        self._check()
        r = { }
        r["selectionpressure"] = self.m_pParams.getSelectionPressure()
        r["elitism"] = self.m_pParams.getUseElitism()
        r["alwaysincludebest"] = self.m_pParams.getAlwaysIncludeBest()
        r["crossoverrate"] = self.m_pParams.getCrossOverRate()

        return r

    def toBytes(self):
        """TODO:"""
        cdef serut.VectorSerializer vSer

        self._check()
        if not self.m_pParams.write(vSer):
            raise InversionParametersException(S(self.m_pParams.getErrorString()))

        return <bytes>vSer.getBufferPointer()[0:vSer.getBufferSize()]

