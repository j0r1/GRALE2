from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cython.operator cimport dereference as deref
import os
import numpy as np
cimport numpy as np
from . import lenses
from . import privutil

include "stringwrappers.pyx"

cimport grale.vector2d as vector2d
cimport grale.multiplanecudacxx as multiplanecudacxx

ctypedef vector2d.Vector2Df Vector2Df

class MultiPlaneCUDAException(Exception):
    """This exception is raised if something goes wrong in the :class:`MultiPlaneCUDA` class."""
    pass

def _getLensParams(lens):
    params = lens.getLensParameters()

    if type(lens) is lenses.CompositeLens:
        returnParams = [ ]
        returnDensity = 0

        for p in params:
            if p["angle"] != 0:
                raise MultiPlaneCUDAException("Currently can't handle rotation in CompositeLens")

            offset = np.array([p["x"], p["y"]], dtype=np.double)
            factor = p["factor"]
            subParams, subDensity = _getLensParams(p["lens"])

            for s in subParams:
                returnParams.append({
                    "position": s["position"] + offset,
                    "width": s["width"],
                    "mass": s["mass"]*factor
                    })

            returnDensity += subDensity*factor

        return returnParams, returnDensity

    elif type(lens) is lenses.PlummerLens:
        return [ { "position": np.zeros((2,), dtype=np.double), "mass": params["mass"], "width": params["width"] }], 0

    elif type(lens) is lenses.MultiplePlummerLens:
        returnParams = [ ]
        for p in params:
            returnParams.append({ "mass": p["mass"], "width": p["width"], 
                                  "position": np.array([p["x"], p["y"]], dtype=np.double) })
        return returnParams, 0

    elif type(lens) is lenses.MassSheetLens:
        density = lens.getLensParameters()["density"]
        return [], density

    raise MultiPlaneCUDAException("Can't handle lens (component) of type {}".format(type(lens)))

defaultLibraryName = "liblens_common.so"

def _findLibraryPath():
    envKey = "GRALE2_MPCUDALIB"
    if envKey in os.environ:
        return os.environ[envKey]

    for p in os.environ["PATH"].split(os.pathsep):
        fullName = os.path.join(p, defaultLibraryName)
        if os.path.exists(fullName):
            return fullName

    raise MultiPlaneCUDAException("Can't find library: environment variable '{}' is not set, and library name '{}' can't be found in PATH".format(envKey, defaultLibraryName))

cdef class MultiPlaneCUDA:
    
    cdef multiplanecudacxx.MultiPlaneCUDA *m_pMPCuda
    cdef double m_angularUnit
    cdef vector[vector[float]] m_massFactors
    cdef vector[float] m_sheetDensities
    cdef cbool m_calculated
    cdef int m_numSources
    cdef object m_plummerParameters
    cdef object m_initialSheetDensities

    def __cinit__(self):
        self.m_pMPCuda = new multiplanecudacxx.MultiPlaneCUDA()
        self.m_angularUnit = 0
        self.m_calculated = False
        self.m_numSources = 0
        self.m_plummerParameters = None

    def __dealloc__(self):
        del self.m_pMPCuda

    def __init__(self, lensesAndRedshifts, thetasAndSourceRedshifts, cosmology="default", libraryPath = None):

        cdef vector[vector[multiplanecudacxx.PlummerInfo]] fixedPlummerParameters
        cdef vector[float] lensRedshifts
        cdef vector[vector[Vector2Df]] rescaledThetas
        cdef vector[float] sourceRedshifts
        cdef int idx, i

        if not libraryPath:
            libraryPath = _findLibraryPath()

        cosmology = privutil.initCosmology(cosmology)
        if not cosmology:
            raise MultiPlaneCUDAException("No cosmological model was specified")

        tmpThetas = [ ]
        for thetas, zs in thetasAndSourceRedshifts:
            sourceRedshifts.push_back(zs)
            thetas = np.array(thetas, dtype=np.double).reshape((-1,2))

            dists = np.sum(thetas**2, 1)**0.5
            distsMax = dists.max()
            if distsMax > self.m_angularUnit:
                self.m_angularUnit = distsMax
    
            tmpThetas.append(thetas)

        self.m_angularUnit /= 10

        # Need to have angular unit here
        idx = 0
        fixedPlummerParameters.resize(len(lensesAndRedshifts))
        self.m_massFactors.resize(len(lensesAndRedshifts))
        self.m_sheetDensities.resize(len(lensesAndRedshifts))

        self.m_plummerParameters = [ ]
        self.m_initialSheetDensities = [ ]
        for lens, zs in lensesAndRedshifts:
            lensRedshifts.push_back(zs)

            params, dens = _getLensParams(lens)
            if len(params) == 0:
                raise MultiPlaneCUDAException("No plummers were specified for lens plane with index {}".format(idx))

            self.m_plummerParameters.append(params)
            self.m_initialSheetDensities.append(dens)

            fixedPlummerParameters[idx].resize(len(params))
            for i in range(len(params)):
                
                par = params[i]
                pos = par["position"]/self.m_angularUnit
                widthInAngularUnit = par["width"]/self.m_angularUnit
                initialMass = par["mass"]
                fixedPlummerParameters[idx][i] = multiplanecudacxx.PlummerInfo(Vector2Df(pos[0], pos[1]),
                                                                               widthInAngularUnit,
                                                                               initialMass)

            self.m_massFactors[idx].resize(fixedPlummerParameters[idx].size())
            idx += 1

        rescaledThetas.resize(len(thetasAndSourceRedshifts))
        for idx in range(len(tmpThetas)):
            thetas = tmpThetas[idx]/self.m_angularUnit

            rescaledThetas[idx].resize(len(thetas))
            for i in range(len(thetas)):
                rescaledThetas[idx][i] = Vector2Df(thetas[i][0], thetas[i][1])

        cp = cosmology.getParameters()
        if not self.m_pMPCuda.init(B(libraryPath), self.m_angularUnit, 
                                   cp["h"], cp["Omega_m"], cp["Omega_r"], cp["Omega_v"], cp["w"],
                                   lensRedshifts, fixedPlummerParameters,
                                   sourceRedshifts, rescaledThetas):
            raise MultiPlaneCUDAException(S(self.m_pMPCuda.getErrorString()))
    
        self.m_numSources = rescaledThetas.size()

    def calculateSourcePositions(self, massFactors, sheetDensities = None):
        cdef int i, j

        if len(massFactors) != self.m_massFactors.size():
            raise MultiPlaneCUDAException("Specified mass factors is for {} lens planes, but expecting {}".format(len(massFactors), self.m_massFactors.size()))

        for i in range(self.m_massFactors.size()):
            if len(massFactors[i]) != self.m_massFactors[i].size():
                raise MultiPlaneCUDAException("Expecting {} factors for lens plane {}, but got {}".format(self.m_massFactors[i].size(), i+1, len(massFactors[i])))
            
            for j in range(self.m_massFactors[i].size()):
                self.m_massFactors[i][j] = massFactors[i][j]

        if sheetDensities is None:
            for i in range(self.m_sheetDensities.size()):
                self.m_sheetDensities[i] = 0;
        else:
            if len(sheetDensities) != self.m_sheetDensities.size():
                raise MultiPlaneCUDAException("Number of mass sheet densities ({}) must equal the number of lens planes ({})".format(len(sheetDensities), self.m_sheetDensities.size()))
            for i in range(self.m_sheetDensities.size()):
                self.m_sheetDensities[i] = sheetDensities[i];

        if not self.m_pMPCuda.calculateSourcePositions(self.m_massFactors, self.m_sheetDensities):
            raise MultiPlaneCUDAException(S(self.m_pMPCuda.getErrorString()))

        self.m_calculated = True

    cdef _sourcePositionConversion(self, const vector[Vector2Df] &positions):
        cdef np.ndarray[double,ndim=2] betas = np.empty([positions.size(),2], dtype = np.double)
        cdef int i

        for i in range(positions.size()):
            betas[i,0] = positions[i].getX()*self.m_angularUnit
            betas[i,1] = positions[i].getY()*self.m_angularUnit
        
        return betas

    cdef _getSourcePositions(self, int srcIdx):
        cdef const vector[Vector2Df] *pPos = self.m_pMPCuda.getSourcePositions(srcIdx)

        if not self.m_calculated:
            raise MultiPlaneCUDAException("Nothing has been calculated yet")
        if pPos == NULL:
            raise MultiPlaneCUDAException(S(self.m_pMPCuda.getErrorString()))
        
        return self._sourcePositionConversion(deref(pPos))

    def getSourcePositions(self, int srcIdx):
        if srcIdx < 0 or srcIdx >= self.m_numSources:
            raise MultiPlaneCUDAException("Invalid source index {} specified, should lie between 0 and {}".format(srcIdx, self.m_numSources-1))

        return self._getSourcePositions(srcIdx)

    def getPlummerParameters(self):
        return self.m_plummerParameters

    def getInitialSheetDensities(self):
        return self.m_initialSheetDensities
