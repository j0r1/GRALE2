from . import constants as CT
from . import inverters
from . import inversionparams
from . import privutil
from . import grid as gridModule
from . import plotutil
import platform
import os
import copy
import random

class InversionException(Exception):
    """TODO:"""
    pass

def estimateStrongLensingMass(Dd, images, skipParamCheck = False):
    """TODO:"""

    count = 0
    massEstimate = 0.0
    
    for imgInfo in images:
        
        img, Ds, Dds = imgInfo["images"], imgInfo["Ds"], imgInfo["Dds"]
        
        if not skipParamCheck:
            params = imgInfo["params"]
            imgType = params["type"]
            if imgType not in [ "extendedimages", "pointimages" ]:
                continue
                
        if img.getNumberOfImages() <= 1:
            continue
        
        tr = img.getTopRightCorner()
        bl = img.getBottomLeftCorner()
        lengthSquared = (((tr-bl)/2)**2).sum()
        
        massEstimate += lengthSquared*((Dd*CT.SPEED_C**2)/(4.0*CT.CONST_G))*(Ds/Dds)
        count += 1
        
    if count == 0:
        raise InversionException("No useful images data sets detected to base mass estimate on")
    
    massEstimate /= count
    return massEstimate        

def _getModuleName(n):
    prefix = "libgalens_"
    suffix = ".so"
    s = platform.system()
    if s == "Windows":
        suffix = ".dll"
        prefix = "galens_"
    elif s == "Darwin":
        suffix = ".dylib"

    return prefix + n + suffix

def _getModuleDirectory(n):
    module = _getModuleName("general")
    modKeyName = "GRALE2_MODULEPATH"
    if modKeyName in os.environ:
        return os.environ[modKeyName]

    condaPrefKey = "CONDA_PREFIX"
    if condaPrefKey in os.environ:
        for s in [ "bin", "lib" ]:
            p = os.path.join(os.environ[condaPrefKey], s)
            n = os.path.join(p, module)
            if os.path.exists(n):
                os.environ[modKeyName] = p # Also set the environment variable
                return p

    # Check if we can find the modules in PATH, or lib dir based on path
    if "PATH" in os.environ:
        subDirs = [ None, [ "..", "lib" ], [ "..", "lib64" ] ]
        for path in os.environ["PATH"].split(os.pathsep):
            for s in subDirs:
                if s is not None:
                    p = os.path.join(path, *s)
                else:
                    p = path
                            
                n = os.path.join(p, module)
                if os.path.exists(n):
                    os.environ[modKeyName] = p # Also set the environment variable
                    return p
            
    raise InversionException("Path with modules not found (GRALE2_MODULEPATH not set, and not detected in other directory)")

def getDefaultModuleParameters(moduleName):
    """TODO:"""

    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set
    return inverters.getInversionModuleDefaultConfigurationParameters(n)

def getDefaultGeneticAlgorithmParameters():
    """TODO:"""
    return inversionparams.GAParameters().getSettings()

def getInversionModuleUsage(moduleName):
    """TODO:"""
    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set
    return inverters.getInversionModuleUsage(n)

def calculateFitness(inputImages, zd, fitnessObjectParameters, lens, moduleName = "general"):
    """TODO:"""
    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set

    # Merge fitnessObjectParameters with defaults
    fullFitnessObjParams = getDefaultModuleParameters(moduleName)
    if fitnessObjectParameters:
        for k in fitnessObjectParameters:
            fullFitnessObjParams[k] = fitnessObjectParameters[k]

    return inverters.calculateFitness(n, inputImages, zd, fullFitnessObjParams, lens)

# TODO: inverters
# TODO: feedback
def invert(inputImages, grid, zd, Dd, popSize, moduleName = "general", massScale = "auto", rescaleBasisFunctions = False, 
           basisFunctionType = "plummer", gridSizeFactor = "default", allowNegativeValues = False, baseLens = None, 
           sheetSearch = "nosheet", fitnessObjectParameters = None, wideSearch = False, maximumGenerations = 16384,
           geneticAlgorithmParameters = { }, inverter = "default", feedbackObject = "default", returnNds = False):
    """TODO:"""

    cellSizeFactorDefaults = { "plummer": 1.7, "gaussian": 1.0, "square": 1.0 }
    print("FeedbackObject1", feedbackObject)
    inverter, feedbackObject = privutil.initInverterAndFeedback(inverter, feedbackObject)
    print("FeedbackObject2", feedbackObject)

    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set
    feedbackObject.onStatus("Full module name is: " + n)
    feedbackObject.onStatus("Detected module path is: " + _getModuleDirectory(n))

    # Merge fitnessObjectParameters with defaults
    fullFitnessObjParams = getDefaultModuleParameters(moduleName)
    if fitnessObjectParameters:
        for k in fitnessObjectParameters:
            fullFitnessObjParams[k] = fitnessObjectParameters[k]

    # Check if we need to convert fractional grid to real grid
    if type(grid) == dict:
        grid = gridModule._fractionalGridToRealGrid(grid)

    # Resize the grid cells if necessary
    if gridSizeFactor == "default":
        cellFactor = cellSizeFactorDefaults[basisFunctionType]
    else:
        cellFactor = gridSizeFactor

    grid = copy.deepcopy(grid) # Make sure we don't modify the input
    for cell in grid:
        cell["size"] *= cellFactor

    # Get massscale
    if massScale == "auto":
        massScale = estimateStrongLensingMass(Dd, inputImages, False)
    elif massScale == "auto_nocheck":
        massScale = estimateStrongLensingMass(Dd, inputImages, True)
    else:
        # Assume the mass scale is a value that needs to be used
        pass

    feedbackObject.onStatus("Mass scale is: {:g} solar masses".format(massScale/CT.MASS_SUN))

    params = inversionparams.GridLensInversionParameters(maximumGenerations, inputImages, grid, Dd, zd, massScale,
               rescaleBasisFunctions, basisFunctionType, allowNegativeValues, baseLens, sheetSearch, fullFitnessObjParams,
               wideSearch)

    # TODO: for now, we're getting the component description in a separate way as it
    #       is not yet integrated in mogal. Once it is, this should be removed
    # By setting the last parameter to None, no calculation is performed and only the
    # fitness components are returned.
    dummyFitness, fitnessComponentDescription = inverters.calculateFitness(n, inputImages, zd, fullFitnessObjParams, None) 

    result = inverter.invert(n, popSize, geneticAlgorithmParameters, params, returnNds)

    # TODO: this should be removed when the TODO above is fixed
    if not returnNds:
        result = (result[0], result[1], fitnessComponentDescription)
    else:
        result = (result[0], fitnessComponentDescription)

    return result

class InversionWorkSpace(object):
    """TODO
    """
    def __init__(self, zLens, cosmology, regionSize, regionCenter = [0, 0], inverter = "default", 
                 renderer = "default", feedbackObject = "default"):
        
        if zLens <= 0 or zLens > 10:
            raise Exception("Invalid lens redshift")
            
        self.zd = zLens
        self.Dd = cosmology.getAngularDiameterDistance(zLens)
        self.imgDataList = []
        self.cosm = cosmology
        self.inversionArgs = { } 

        if regionSize <= 0:
            raise Exception("Invalid region size")
        
        # Rough indication of grid position and dimensions, but randomness
        # may be added unless specified otherwise
        self.regionSize = regionSize
        self.regionCenter = copy.deepcopy(regionCenter)
        self.grid = None
        
        self.renderer = renderer
        self.inverter = inverter
        self.feedbackObject = feedbackObject

    def getCosmology(self):
        return self.cosm
        
    def clearImageDataList(self):
        self.imgDataList = []

    def getImageDataList(self):
        return self.imgDataList
        
    def addImageDataToList(self, imgDat, zs, imgType, otherParameters = {}):
        # check that imgDat exists
        num = imgDat.getNumberOfImages()
        
        if zs < self.zd:
            raise Exception("Can't add a source with a smaller redshift than the lens")
        
        params = copy.deepcopy(otherParameters)
        if imgType:
            params["type"] = imgType
            
        entry = {
            "images": imgDat,
            "Ds": self.cosm.getAngularDiameterDistance(zs),
            "Dds": self.cosm.getAngularDiameterDistance(self.zd, zs),
            "params": params
        }
        self.imgDataList.append(entry)
        
    # Overrides the grid
    def setGrid(self, grid):
        self.grid = grid

    def getGrid(self):
        return self.grid

    def _getGridDimensions(self, randomFraction):
        # TODO: make it possible to specify your own function for the randomness?
        w = self.regionSize
        dx = (random.random()-0.5) * w*randomFraction
        dy = (random.random()-0.5) * w*randomFraction
        c = [ self.regionCenter[0], self.regionCenter[1]]
        return w, c
    
    def setUniformGrid(self, subDiv, randomFraction = 0.05):
        w, c = self._getGridDimensions(randomFraction)
        self.grid = gridModule.createUniformGrid(w, c, subDiv)

        # TODO: onstatus stuff?
        # TODO: callback to intercept new grid?
        
    def setSubdivisionGrid(self, lensOrLensInfo, minSquares, maxSquares, startSubDiv = 1, randomFraction = 0.05):
        w, c = self._getGridDimensions(randomFraction)
        if type(lensOrLensInfo) == dict:
            lensInfo = lensOrLensInfo
        else:
            extraFrac = 1.0001
            lensInfo = { 
                "lens": lensOrLensInfo,
                # Make it slightly wider to avoid going out of bounds
                "bottomleft": [ c[0] - (w*extraFrac)/2, c[1] - (w*extraFrac)/2 ],
                "topright": [ c[0] + (w*extraFrac)/2, c[1] + (w*extraFrac)/2 ],
            }
            
        # Make sure we have the density points calculated, abuse one of the plot functions
        # for this
        plotutil.plotDensity(lensInfo, axes = False)
        
        self.grid = gridModule.createSubdivisionGrid(w, c, lensInfo, minSquares, maxSquares, startSubDiv)
       
    def setDefaultInversionArguments(self, **kwargs):
        self.inversionArgs = kwargs

    def invert(self, populationSize, **kwargs):
        newKwargs = { }
        newKwargs["inverter"] = self.inverter
        newKwargs["feedbackObject"] = self.feedbackObject
        for a in self.inversionArgs:
            newKwargs[a] = self.inversionArgs[a]

        for a in kwargs:
            newKwargs[a] = kwargs[a]

        lens = invert(self.imgDataList, self.grid, self.zd, self.Dd, populationSize, **newKwargs)
        return lens

    def calculateFitness(self, lens):
        fitnessObjectParameters = None
        moduleName = "general"
        if "fitnessObjectParameters" in self.inversionArgs:
            fitnessObjectParameters = self.inversionArgs["fitnessObjectParameters"]
        if "moduleName" in self.inversionArgs:
            moduleName = self.inversionArgs["moduleName"]

        return calculateFitness(self.imgDataList, self.zd, fitnessObjectParameters, lens, moduleName)

def getDefaultInverter():
    return inverters.getDefaultInverter()

def setDefaultInverter(x):
    inverters.setDefaultInverter(x)

