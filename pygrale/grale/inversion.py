from . import constants as CT
from . import inverters
from . import inversionparams
from . import privutil
from .grid import _fractionalGridToRealGrid
import platform
import os
import copy

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

# TODO: inverters
# TODO: feedback
def invert(inputImages, grid, zd, Dd, popSize, moduleName = "general", massScale = "auto", rescaleBasisFunctions = False, 
           basisFunctionType = "plummer", gridSizeFactor = "default", allowNegativeValues = False, baseLens = None, 
           sheetSearch = "nosheet", fitnessObjectParameters = None, wideSearch = False, maximumGenerations = 16384,
           geneticAlgorithmParameters = { }, inverter = "singlecore", feedbackObject = "default"):
    """TODO:"""

    cellSizeFactorDefaults = { "plummer": 1.7, "gaussian": 1.0, "square": 1.0 }
    inverter, feedbackObject = privutil.initInverterAndFeedback(inverter, feedbackObject)

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
        grid = _fractionalGridToRealGrid(grid)

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

    resultLens = inverter.invert(n, popSize, geneticAlgorithmParameters, params)

    return resultLens
