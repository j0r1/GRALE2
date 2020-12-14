"""This module contains functions for gravitational lens inversion. The
most straightforward method is to use the class :class:`InversionWorkSpace`,
which keeps track of the settings and provides an easier interface to
the other functions. For absolute control however, the individual functions
can still be used.

This module will try do determine automatically where the inversion libraries
to be used in the genetic algorithm are to be found. If the environment
variable ``GRALE2_MODULEPATH`` is set, then this path will always be used.
"""
from __future__ import print_function
from . import constants as CT
from . import inverters
from . import inversionparams
from . import privutil
from . import grid as gridModule
from . import plotutil
from . import lenses
from . import multiplane
import platform
import os
import copy
import random
import numpy as np

_cellSizeFactorDefaults = { "plummer": 1.7, "gaussian": 1.0, "square": 1.0 }

class InversionException(Exception):
    """An exception that's generated if something goes wrong in a function
    provided by this module."""
    pass

def estimateStrongLensingMass(Dd, images, skipParamCheck = False):
    """For a gravitational lens that's located at angular diameter distance
    `Dd`, estimate the strong lensing mass based on the images data in
    `images`. Each entry in the `images` list should be a dictionary with
    at least these keys:

     - ``images``:  the :class:`ImagesData<grale.images.ImagesData>` instance
       describing the images of a source
     - ``Ds`` and ``Dds``: angular diameter distances to this source and from
       lens to source respectively.

    Unless `skipParamCheck` is ``True``, the dictionary should also have a 
    ``params`` key, for futher parameters of this images data set. These
    ``params`` should also be a dictionary with at least a ``type`` field;
    the images data set is skipped unless this type field is either
    ``extendedimages``, ``pointimages`` or ``pointgroupimages``. This way, images data sets describing
    null space information for example, are skipped.

    If there's only one image, it will not be used in the mass estimate. Otherwise,
    the mass estimate for a single images data set will be estimated as the point
    mass that would be needed to create an Einstein ring with the same diameter as
    the images' separaration. The final mass estimate is the average of these
    separate estimates.
    """

    count = 0
    massEstimate = 0.0
    
    for imgInfo in images:
        
        img, Ds, Dds = imgInfo["images"], imgInfo["Ds"], imgInfo["Dds"]
        
        if not skipParamCheck:
            params = imgInfo["params"]
            imgType = params["type"]
            if imgType not in [ "extendedimages", "pointimages", "pointgroupimages" ]:
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

def getDefaultModuleParameters(moduleName = "general"):
    """For the specified module name, query the default parameters."""

    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set
    return inverters.getInversionModuleDefaultConfigurationParameters(n)

def getDefaultGeneticAlgorithmParameters():
    """Returns the default parameters for the genetic algorithm that's used in lens inversion."""
    return inversionparams.GAParameters().getSettings()

def getInversionModuleUsage(moduleName = "general"):
    """Returns a usage description that's provided by the specified genetic algorithm 
    module. The usage information for the ``general`` module can be viewed here: `usage <./usage_general.html>`_
    """
    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set
    return inverters.getInversionModuleUsage(n)

def _mergeModuleParameters(fitnessObjectParameters, moduleName, cosmology):
    fullFitnessObjParams = getDefaultModuleParameters(moduleName)
    if fitnessObjectParameters:
        for k in fitnessObjectParameters:
            fullFitnessObjParams[k] = fitnessObjectParameters[k]

    if "general_cosmology" in fullFitnessObjParams:
        if cosmology is None:
            fullFitnessObjParams["general_cosmology"] = None
        else:
            p = cosmology.getParameters()
            fullFitnessObjParams["general_cosmology"] = np.array([
                p["h"], p["Omega_m"], p["Omega_r"], p["Omega_v"], p["w"]
            ])
    
    return fullFitnessObjParams

def calculateFitness(inputImages, zd, fitnessObjectParameters, lensOrBackProjectedImages, moduleName = "general", cosmology = None):
    """You can pass several parameters in the same way as you would do for a
    lens inversion, but here you also specify the specific lens for which the
    relevant fitness measures should be calculated or the pre-calculated
    back-projected images.

    Arguments:
     - `inputImages`: list of input images data instances that should be used in the
       inversion. See the :func:`invert` function for more information about the
       format.
     - `zd`: the redshift to the lens.
     - `fitnessObjectParameters`: parameters to be used when initializing the inversion
       module that's specified
     - `lensOrBackProjectedImages`: either the gravitational lens model for which the 
       fitness measures should be calculated, or a similar list as `inputImages`,
       for which the positions have already been mapped onto the source plane.
     - `moduleName`: name of the inversion module for the genetic algorithm.
    """
    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set

    # Merge fitnessObjectParameters with defaults
    fullFitnessObjParams = _mergeModuleParameters(fitnessObjectParameters, moduleName, cosmology)

    lens, bpImages = None, None
    try:
        d = lensOrBackProjectedImages.getLensDistance()
        lens = lensOrBackProjectedImages
    except:
        pass

    try:
        l = len(lensOrBackProjectedImages)
        bpImages = lensOrBackProjectedImages
    except:
        pass

    return inverters.calculateFitness(n, inputImages, zd, fullFitnessObjParams, lens=lens, bpImages=bpImages)

def _invertCommon(inverter, feedbackObject, moduleName, fitnessObjectParameters,
                  massScale, DdAndZd, inputImages, getParamsFunction,
                  popSize, geneticAlgorithmParameters, returnNds, cosmology):

    # Set to some bad values as they don't make sense for a multi-plane inversion
    Dd, zd = DdAndZd if DdAndZd else (None, float("NaN"))

    inverter, feedbackObject = privutil.initInverterAndFeedback(inverter, feedbackObject)

    n = _getModuleName(moduleName)
    _getModuleDirectory(n) # So that GRALE2_MODULEPATH gets set
    feedbackObject.onStatus("Full module name is: " + n)
    feedbackObject.onStatus("Detected module path is: " + _getModuleDirectory(n))

    # Merge fitnessObjectParameters with defaults
    fullFitnessObjParams = _mergeModuleParameters(fitnessObjectParameters, moduleName, cosmology)

    # Get massscale
    if massScale == "auto":
        massScale = estimateStrongLensingMass(Dd, inputImages, False)
    elif massScale == "auto_nocheck":
        massScale = estimateStrongLensingMass(Dd, inputImages, True)
    else:
        # Assume the mass scale is a value that needs to be used
        pass

    feedbackObject.onStatus("Mass scale is: {:g} solar masses".format(massScale/CT.MASS_SUN))

    params = getParamsFunction(fullFitnessObjParams, massScale)

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

def invertMultiPlane(cosmology, inputImages, basisLensesAndRedshifts, popSize, moduleName="general_gpu",
                     massScale="auto", allowNegativeValues=False, sheetSearch="nosheet",
                     fitnessObjectParameters=None, massScaleSearchType="regular",
                     maximumGenerations=16384, geneticAlgorithmParameters={ },
                     returnNds=False, deviceIndex = "rotate",
                     inverter="default", feedbackObject="default"):
    """Perform a multi-plane lens inversion. This is a rather low-level function,
    it may be easier to use an instance of :class:`InversionWorkSpace` instead.
    
    Arguments:

     - `cosmology`: The cosmological model to use, to calculate the necessary angular
       diameter distances based on the specified redshifts.

     - `inputImages`: a list of dictionaries with the following entries:
     
        - ``images``: an :class:`ImagesData<grale.images.ImagesData>` instance that
          describes the images of a source, the null space etc.
        - ``z``: the redshift of this source.
        - ``params``: for the ``general_gpu`` module, this could contain a 
          dictionary with at least a ``type`` field, which can be e.g. 
          ``pointimages`` or ``extendednullgrid`` (see the `usage <./usage_general.html>`_ 
          documentation). 

     - `basisLensesAndRedshifts`: A list representing the different lens planes.
       Each list entry is a dictionary with needs a key `z` to describe the
       plane's redshift, and a key `lenses` that describes the basis functions
       in this lens plane. This `lenses` entry should be a list of dictionaries
       with the following entries:

         - `lens`: the lens model for this basis function.
         - `center`: the x,y position at which this lens model should be placed.
         - `mass`: the mass of this lens model, in the relevant area. Some models have a total
           mass parameter which would likely work fine, but not all models have this (e.g. a
           SIS lens or a mass sheet). You should then precalculate the mass in the strong lensing
           region (approximately) and store it in this entry. This is needed to the algorithm can
           estimate the total lensing mass of a certain combination of weighted basis functions.

     - `popSize`: the size of the population in the genetic algorithm, e.g. 512.

     - `moduleName`: name of the inversion module for the genetic algorithm.

     - `massScale`: a rough estimate of the total mass of the gravitational lens.
       Set to ``"auto"`` or ``"auto_nocheck"`` to let the :func:`estimateStrongLensingMass`
       function provide this estimate automatically.

     - `allowNegativeValues`: by default, the weight of the basis functions are only
       allowed to be positive, to make certain that an overall positive mass density
       is obtained. In case corrections to a certain mass distribution are being sought,
       negative weights can be allowed by setting this parameter to ``False``.

     - `sheetSearch`: by default, only the basis functions for the grid cells are used.
       You can also allow a mass-sheet basis function, which may be useful as this kind
       of effect is difficult to model by a grid of basis functions. To do so, set this
       to ``"genome"``, in which a mass-sheet component will be added to each lens plane.

     - `fitnessObjectParameters`: parameters for the lens inversion module for the
       generic algorithm. For the ``"general"`` module, more information can be
       found in the `usage <./usage_general.html>`_ documentation.

     - `massScaleSearchType`: by default (``"regular"``), a relatively narrow mass 
       range around the provided mass estimate will be explored. To make this search 
       wider (which can be useful if you're less certain of the total mass, e.g. 
       when including weak lensing measurements over a larger area), you can set this
       parameter to ``"wide"``. It can also be set to ``"nosearch"`` to disable the
       mass scale search completely. Finally, mainly for testing purposes, it can also
       be set to a dictionary with the following entries:

        - `startFactor`
        - `stopFactor`
        - `numIterations`
        - `firstIterationSteps`
        - `nextIterationSteps`

     - `maximumGenerations`: if the genetic algorithm didn't stop by itself after
       this many generations, stop it anyway. To test an inversion script completely,
       it can be useful to temporarily stop the genetic algorithm after only a small
       number of generations so that the code doesn't take long to run.

     - `geneticAlgorithmParameters`: a dictionary with general genetic algorithm parameters
       that should be changed from their defaults. Known names and their defaults are

        - ``selectionpressure`` (default is 2.5)
        - ``elitism`` (default is ``True``)
        - ``alwaysincludebest`` (default is ``True``)
        - ``crossoverrate`` (default is 0.9)

       For more information about their meaning, refer to the `documentation <http://research.edm.uhasselt.be/jori/mogal/documentation/classmogal_1_1GeneticAlgorithmParams.html>`_
       of the library that's used for the genetic algorithm.

     - `returnNds`: by default, this function will return a single gravitational lens
       model. If there are several fitness measures however, the end result is actually
       a non-dominated set of models. The inversion module for the genetic algorithm has
       some default strategy for choosing one solution from this set. In case you'd like
       to get the complete non-dominated set instead, you can set this flag to ``True``.

     - `deviceIndex`: this multi-plane inversion uses a GPU to back-project the image
       data, and by setting a specific number, a specific device can be specified. To
       allow multiple GPUs to be used automatically, you can leave this to ``"auto"``
       and use an :mod:`inverter <grale.inverters>` with as many processes as you have
       GPUs.

     - `inverter`: specifies the inverter to be used. See the :mod:`inverters<grale.inverters>`
       module for more information.

     - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.
    """
        
    if massScale == "auto" or massScale == "auto_nocheck":
        minZd = min([entry["z"] for entry in basisLensesAndRedshifts])
        minDd = cosmology.getAngularDiameterDistance(minZd)

        # We're going to use a list with Dds entries based on the first lens plane,
        # for mass determination
        newImages = []
        for entry in inputImages:
            entry = entry.copy()
            entry["Ds"] = cosmology.getAngularDiameterDistance(entry["z"])
            entry["Dds"] = cosmology.getAngularDiameterDistance(minZd, entry["z"])
            newImages.append(entry)

        massScale = estimateStrongLensingMass(minDd, newImages, False if massScale == "auto" else True)

    def getParamsFunction(fullFitnessObjParams, massScale):
        return inversionparams.LensInversionParametersMultiPlaneGPU(cosmology,
                basisLensesAndRedshifts, inputImages, massScale, sheetSearch,
                fullFitnessObjParams, maximumGenerations, allowNegativeValues,
                massScaleSearchType, deviceIndex)

    return _invertCommon(inverter, feedbackObject, moduleName, fitnessObjectParameters,
                  massScale, None, inputImages, getParamsFunction, popSize,
                  geneticAlgorithmParameters, returnNds, cosmology)

def invert(inputImages, basisFunctions, zd, Dd, popSize, moduleName = "general", massScale = "auto",
           allowNegativeValues = False, baseLens = None, 
           sheetSearch = "nosheet", fitnessObjectParameters = None, massScaleSearchType = "regular", maximumGenerations = 16384,
           geneticAlgorithmParameters = { }, returnNds = False, inverter = "default", feedbackObject = "default",
           cosmology = None):
    """Start the genetic algorithm to look for a gravitational lens model that's
    compatible with the specified input images. This is a rather low-level function,
    it may be easier to use an instance of :class:`InversionWorkSpace` instead.
    This function is for a single lens plane inversion, :func:`invertMultiPlane` is
    the multi-lens plane counterpart.

    Arguments:
     - `inputImages`: a list of dictionaries with the following entries:
     
        - ``images``: an :class:`ImagesData<grale.images.ImagesData>` instance that
          describes the images of a source, the null space etc.
        - ``Ds`` and ``Dds``: the angular diameter distances to this source.
        - ``params``: not used for older inversion modules for the genetic
          algorithm, but for the ``general`` module, this should contain a 
          dictionary with at least a ``type`` field, which can be e.g. 
          ``pointimages`` or ``extendednullgrid`` (see the `usage <./usage_general.html>`_ 
          documentation). Other parameters may be set as well, e.g. you could set
          ``timedelay`` to ``False`` to ignore the time delay information in
          a particular images data instance.

     - `basisFunctions`: a list of entries containing basis function info.
       This parameter should be a list of dictionaries with the following entries:

         - `lens`: the lens model for this basis function.
         - `center`: the x,y position at which this lens model should be placed.
         - `mass`: the mass of this lens model, in the relevant area. Some models have a total
           mass parameter which would likely work fine, but not all models have this (e.g. a
           SIS lens or a mass sheet). You should then precalculate the mass in the strong lensing
           region (approximately) and store it in this entry. This is needed to the algorithm can
           estimate the total lensing mass of a certain combination of weighted basis functions.

     - `zd`: the redshift to the lens, used when calculating time delays
       (for other purposes the angular diameter distance will be used).

     - `Dd`: the angular diameter distance to the lens.

     - `popSize`: the size of the population in the genetic algorithm, e.g. 512.

     - `moduleName`: name of the inversion module for the genetic algorithm.

     - `massScale`: a rough estimate of the total mass of the gravitational lens.
       Set to ``"auto"`` or ``"auto_nocheck"`` to let the :func:`estimateStrongLensingMass`
       function provide this estimate automatically.

     - `allowNegativeValues`: by default, the weight of the basis functions are only
       allowed to be positive, to make certain that an overall positive mass density
       is obtained. In case corrections to a certain mass distribution are being sought,
       negative weights can be allowed by setting this parameter to ``False``.

     - `baseLens`: in case an approximate solution is already known, it can be set as
       the lens to start from. In this scenario, the `allowNegativeValues` parameter
       is usually set to ``True``. Note that when the base lens is specified, only the
       modifications to it are returned by the inversion routine, and if you want the
       full lens you need to combine this base lens and the modifications by using a
       :class:`CompositeLens <grale.lenses.CompositeLens>` instance.

     - `sheetSearch`: by default, only the basis functions for the grid cells are used.
       You can also allow a mass-sheet basis function, which may be useful as this kind
       of effect is difficult to model by a grid of basis functions. The traditional
       way to specify that a mass sheet basis function should be used is to set this
       to ``"genome"``. Alternatively you can specify a :mod:`lens model <grale.lenses>`
       to be used in a similar way, allowing you to e.g. use a mass disk. Of course,
       you can use any lens model here, it doesn't need to be one resembling a mass
       sheet. The lens model is very roughly an upper limit in the sense that the
       search for this weight is initially between 0 and 1 (it is not included in the
       scaling loop used in the genetic algorithm).

     - `fitnessObjectParameters`: parameters for the lens inversion module for the
       generic algorithm. For the ``"general"`` module, more information can be
       found in the `usage <./usage_general.html>`_ documentation.

     - `massScaleSearchType`: by default (``"regular"``), a relatively narrow mass 
       range around the provided mass estimate will be explored. To make this search 
       wider (which can be useful if you're less certain of the total mass, e.g. 
       when including weak lensing measurements over a larger area), you can set this
       parameter to ``"wide"``. It can also be set to ``"nosearch"`` to disable the
       mass scale search completely. Finally, mainly for testing purposes, it can also
       be set to a dictionary with the following entries:

        - `startFactor`
        - `stopFactor`
        - `numIterations`
        - `firstIterationSteps`
        - `nextIterationSteps`

     - `maximumGenerations`: if the genetic algorithm didn't stop by itself after
       this many generations, stop it anyway. To test an inversion script completely,
       it can be useful to temporarily stop the genetic algorithm after only a small
       number of generations so that the code doesn't take long to run.

     - `geneticAlgorithmParameters`: a dictionary with general genetic algorithm parameters
       that should be changed from their defaults. Known names and their defaults are

        - ``selectionpressure`` (default is 2.5)
        - ``elitism`` (default is ``True``)
        - ``alwaysincludebest`` (default is ``True``)
        - ``crossoverrate`` (default is 0.9)

       For more information about their meaning, refer to the `documentation <http://research.edm.uhasselt.be/jori/mogal/documentation/classmogal_1_1GeneticAlgorithmParams.html>`_
       of the library that's used for the genetic algorithm.

     - `returnNds`: by default, this function will return a single gravitational lens
       model. If there are several fitness measures however, the end result is actually
       a non-dominated set of models. The inversion module for the genetic algorithm has
       some default strategy for choosing one solution from this set. In case you'd like
       to get the complete non-dominated set instead, you can set this flag to ``True``.

     - `inverter`: specifies the inverter to be used. See the :mod:`inverters<grale.inverters>`
       module for more information.

     - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

     - `cosmology`: depending on the inversion algorithm, it may be necessary to specify a
       cosmological model.
    """

    def getParamsFunction(fullFitnessObjParams, massScale):
        return inversionparams.LensInversionParametersSinglePlaneCPU(maximumGenerations, inputImages, basisFunctions,
                                                         Dd, zd, massScale, allowNegativeValues, baseLens, 
                                                         sheetSearch, fullFitnessObjParams, massScaleSearchType)

    return _invertCommon(inverter, feedbackObject, moduleName, fitnessObjectParameters,
                  massScale, [Dd, zd], inputImages, getParamsFunction, popSize,
                  geneticAlgorithmParameters, returnNds, cosmology)

def defaultLensModelFunction(operation, operationInfo, parameters):
    """This is the default `lensModelFunction` that's used in 
    :func:`InversionWorkSpace.addBasisFunctionsBasedOnCurrentGrid <grale.inversion.InversionWorkSpace.addBasisFunctionsBasedOnCurrentGrid>`
    This default function will case the same behaviour as the 
    classic grid based inversion procedure.

    The `initialParameters` in that function call will be merged with
    a dictionary with the following entries:
    
     - ``basistype``: defaults to ``plummer``, but can also be ``square`` or ``gaussian``
     - ``sizefactor``: a specific value can be set here which converts the grid square
       size to a width of the basis function. The default value depends on the basis
       function type: 1.7 ``for plummer`` and 1.0 for the two others.
     - ``rescale``: defaults to ``False``, and has the same meaning as
       `rescaleBasisFunctions` in the documentation of :func:`invert`.
     - ``totalmass``: the mass scale for the entire lensing region. In case it's
       ``"auto"``, it will be estimated from the images stored in the
       :class:`InversionWorkSpace` instance.
    """

    if operation == "start":

        # Initialize to default parameters
        usedParams = { "basistype": "plummer", "sizefactor": "default", "rescale": False, 
                       "totalmass": "auto" }
        if parameters:
            for k in parameters:
                usedParams[k] = parameters[k]

        iws = operationInfo["workspace"]
        lpIdx = operationInfo["lensplaneindex"]
        totalMass = iws.estimateStrongLensingMass() if usedParams["totalmass"] == "auto" else usedParams["totalmass"]

        if usedParams["sizefactor"] == "default":
            usedParams["sizefactor"] = _cellSizeFactorDefaults[usedParams["basistype"]]

        grid = gridModule._fractionalGridToRealGrid(operationInfo["grid"])
        if not usedParams["rescale"]:
            numCells = len(grid)
            usedParams["cellmass"] = totalMass/numCells
        else:
            totalRescale = 0
            for cell in grid:
                totalRescale += (cell["size"]*usedParams["sizefactor"]/CT.ANGLE_ARCSEC)**2
            usedParams["cellmass"] = totalMass/totalRescale

        usedParams["Dd"] = iws.getLensDistance(lpIdx)
        return usedParams

    if operation != "add":
        raise InversionException(f"Unexpected operationtype {operation}, was expecting 'add'")

    width = operationInfo["size"]*parameters["sizefactor"]
    mass = parameters["cellmass"]
    if parameters["rescale"]:
        mass *= (width/CT.ANGLE_ARCSEC)**2 # smaller cells have less mass
    
    classes = { "plummer": lenses.PlummerLens,
                "gaussian": lenses.GaussLens,
                "square": lenses.SquareLens }

    lens = classes[parameters["basistype"]](parameters["Dd"], { "width": width, "mass": mass })
    return lens, mass


class InversionWorkSpace(object):
    """This class tries to make it more straightforward to perform lens inversions
    by automatically keeping track of the inversion region, the input images and
    several other settings.

    Typically, you'd perform the following steps:

     - create an instance of the class, specifying the redshift to the lens, the
       cosmological model and the region for the inversion.

     - add the images data sets to be used in the inversion using calls to
       :func:`addImageDataToList`.

     - create a uniform grid for a first inversion with :func:`setUniformGrid` and
       run the inversion with a call to :func:`invert`.

     - based on the returned lens model, create a grid with smaller cells in
       regions with more mass using :func:`setSubdivisionGrid`, and again run the
       inversion with a call to :func:`invert`.

     - repeat as needed/desired.

    To get more control over the basis functions of which the weights will be
    optimized, you can also set basis functions directly (see e.g. :func:`addBasisFunctions`)
    and optimize their weights using :func:`invertBasisFunctions`.
    """
    def __init__(self, zLens, regionSize, regionCenter = [0, 0], inverter = "default", 
                 renderer = "default", feedbackObject = "default", cosmology = "default"):
        
        """Constructor for this class.

        Arguments:
         - `zLens`: the redshift to the gravitational lens for a single lens plane inversion,
           or a list of redshifts for an (experimental) multi lens plane inversion.

         - `regionSize` and `regionCenter`: the width and height of the region in which the
           inversion should take plane, as well as its center. This will be used to base the
           grid dimensions on, but by default some randomness will be added (see e.g. :func:`setUniformGrid`).
           
           These are the default values in :func:`setUniformGrid` and :func:`setSubdivisionGrid`
           but can still be overridden there.

           In case of a multi-plane inversion, you can specify different sizes for each plane
           if desired, otherwise the same settings will be used for all planes.

         - `inverter`: specifies the inverter to be used. See the :mod:`inverters<grale.inverters>`

         - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
           to speed up the calculation of the mass densities (used for the procedure with the
           subdivision grid)
        
         - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

         - `cosmology`: an instance of :class:`Cosmology <grale.cosmology.Cosmology>`, describing
           the cosmological model that should be used throughout these inversions.
        """

        if type(zLens) == list:
            prevZ = zLens[0]
            for z in zLens[1:]:
                if z <= prevZ:
                    raise InversionException("Lens plane redshifts need to be strictly increasing")
                prevZ = z

            self.isMultiPlane = True
        else:
            self.isMultiPlane = False
            zLens = [ zLens ]

        for z in zLens:
            if z <= 0 or z > 10:
                raise InversionException("Invalid lens redshift {}".format(z))

        cosmology = privutil.initCosmology(cosmology)
        if not cosmology:
            raise InversionException("No cosmological model was specified")
            
        self.zd = zLens
        self.Dd = [ cosmology.getAngularDiameterDistance(z) for z in zLens ]
        self.imgDataList = []
        self.cosm = cosmology
        self.inversionArgs = { }

        if type(regionSize) != list:
            regionSize = [ copy.deepcopy(regionSize) for x in range(len(zLens)) ]
            regionCenter = [ copy.deepcopy(regionCenter) for x in range(len(zLens)) ]

        if len(regionSize) != len(zLens):
            raise InversionException("There must be the same amount of region sizes as lens planes")

        for r in regionSize:
            if r <= 0:
                raise InversionException("Invalid region size")
        
        # Rough indication of grid position and dimensions, but randomness
        # may be added unless specified otherwise
        self.regionSize = copy.deepcopy(regionSize)
        self.regionCenter = copy.deepcopy(regionCenter)
        self.grid = [ None for z in zLens ]
        self.clearBasisFunctions() # initializes self.basisFunctions
        
        self.renderer = renderer
        self.inverter = inverter
        self.feedbackObject = feedbackObject

    def _checkLensPlaneIndex(self, lpIdx):
        if self.isMultiPlane:
            if lpIdx is None:
                raise InversionException("For a multi-plane scenario, a lens plane index is required")
            if lpIdx < 0 or lpIdx >= len(self.zd):
                raise InversionException("Invalid index {} for lens plane".format(lpIdx))

            return lpIdx

        # Single plane
        if lpIdx is None:
            return 0
        if lpIdx != 0:
            raise InversionException("Only lens plane index 0 is allowed in single lens plane scenario")
        return lpIdx

    def getLensDistance(self, lpIdx = None):
        """Returns the angular diameter distance to the lens. In case a multi-plane
        setting is used, ``lpIdx`` needs to specify a particular lens plane."""
        lpIdx = self._checkLensPlaneIndex(lpIdx)
        return self.Dd[lpIdx]

    def setRegionSize(self, regionSize, regionCenter = [0, 0], lpIdx = "all"):
        """Set the inversion region, same as in the constructor: `regionSize` and 
        `regionCenter` are the width, height and center of the region in which the
        inversion should take plane. This will be used to base the
        grid dimensions on, but by default some randomness will be added 
        (see e.g. :func:`setUniformGrid`).
        
        In a multi-plane scenario, ``lpIdx`` can be used to specify a particular
        lens plane for which the region should be set.
        """
        if lpIdx == "all":
            for i in range(len(self.regionSize)):
                self.regionSize[i] = copy.deepcopy(regionSize)
                self.regionCenter[i] = copy.deepcopy(regionCenter)
        else:
            lpIdx = self._checkLensPlaneIndex(lpIdx)
            self.regionSize[lpIdx] = regionSize
            self.regionCenter[lpIdx] = copy.deepcopy(regionCenter)

    def getCosmology(self):
        """Returns the cosmological model that was specified during initialization."""
        return self.cosm
        
    def clearImageDataList(self):
        """Clears the list of images data information."""
        self.imgDataList = []

    def getImageDataList(self):
        """Returns the list that contains the images data sets, and which is built up
        by successive calls to :func:`addImageDataToList`."""
        return self.imgDataList
        
    def addImageDataToList(self, imgDat, zs, imgType, otherParameters = {}):
        """Adds the :class:`ImagesData <grale.images.ImagesData>` instance in `imgDat`
        to the list of inputs. The corresponding redshift is `zs`, and `imgType`
        describes the type of input, e.g. ``"extendedimages"``. The
        `otherParameters` dictionary provides additional settings for this
        input set. See the `usage <./usage_general.html>`_ documentation for other types
        and parameters.
        """
        # check that imgDat exists
        num = imgDat.getNumberOfImages()
        
        if zs is None: # TODO: document this
            Dds = 1.0
            Ds = 1.0
        else:
            if zs < self.zd[0]:
                raise InversionException("Can't add a source with a smaller redshift than the closest lens")
            Ds = self.cosm.getAngularDiameterDistance(zs)
            # TODO: for backward compatibility, we'll use the first lens plane here
            Dds = self.cosm.getAngularDiameterDistance(self.zd[0], zs)
        
        params = copy.deepcopy(otherParameters)
        if imgType:
            params["type"] = imgType
            
        # TODO: clean this so that at this point only the redshifts are stored
        entry = {
            "images": imgDat,
            "Ds": Ds,
            "Dds": Dds,
            "z": zs,
            "params": params
        }
        self.imgDataList.append(entry)
        
    # Overrides the grid
    def setGrid(self, grid, lpIdx = None):
        """Usually, the functions :func:`setUniformGrid` and :func:`setSubdivisionGrid`
        will be used to control the grid (which in turn controls the layout of the
        basis functions). If this doesn't suffice, you can provide a specific grid
        obtained by one of the functions in :mod:`grid <grale.grid>` yourself using 
        this function.
        
        In a multi-plane setting, ``lpIdx`` specifies the lens plane for which this
        grid is relevant.
        """
        lpIdx = self._checkLensPlaneIndex(lpIdx)
        self.grid[lpIdx] = grid

    def getGrid(self, lpIdx = None):
        """Retrieves the currently set grid, e.g. for plotting using 
        :func:`plotSubdivisionGrid <grale.plotutil.plotSubdivisionGrid>`.
        
        In a multi-plane setting, ``lpIdx`` specifies a particular lens
        plane.
        """
        lpIdx = self._checkLensPlaneIndex(lpIdx)
        return self.grid[lpIdx]

    def _getGridDimensions(self, randomFraction, regionSize, regionCenter, lpIdx):
        if regionSize is None:
            regionSize = self.regionSize[lpIdx]
        if regionCenter is None:
            regionCenter = self.regionCenter[lpIdx]

        if type(randomFraction) == float: # Just a number, use the default method
            w = regionSize
            dx = (random.random()-0.5) * w*randomFraction
            dy = (random.random()-0.5) * w*randomFraction
            c = [ regionCenter[0] + dx, regionCenter[1] + dy ]
        else: # assume it's something we can call to obtain the grid dimensions
            w, c = randomFraction(regionSize, copy.copy(regionCenter))

        return w, c
    
    def setUniformGrid(self, subDiv, randomFraction = 0.05, regionSize = None, regionCenter = None,
                       lpIdx = "all"):
        """Based on the size and center provided during initialization, create
        a uniform grid with `subDiv` subdivisions along each axis, resulting
        in `subDiv`x`subDiv` cells. 
        
        The `randomFraction` parameter controls some randomness in the grid
        positioning: if this is a number, then the center will be offset randomly
        in X and Y directions by a fraction of the grid size. You can also specify
        a function instead of a number. In that case, it will receive the region
        size and center as two arguments, and it should return a tuple
        containing the width and center of the actual grid to use. This gives
        you somewhat more freedom in adding randomness to the grid size and center.

        If specified, `regionSize` and `regionCenter` override the internally stored
        dimensions.

        In a multi-plane scenario, ``lpIdx`` can be used to specify only one
        lens plane.
        """
        
        if lpIdx == "all":
            for i in range(len(self.grid)):
                w, c = self._getGridDimensions(randomFraction, regionSize, regionCenter, i)
                self.grid[i] = gridModule.createUniformGrid(w, c, subDiv)
        else:
            lpIdx = self._checkLensPlaneIndex(lpIdx)
            w, c = self._getGridDimensions(randomFraction, regionSize, regionCenter, lpIdx)
            self.grid[lpIdx] = gridModule.createUniformGrid(w, c, subDiv)

    def _getSinglePlaneSubDivGrid(self, lensOrLensInfo,  minSquares, maxSquares, startSubDiv, randomFraction, regionSize, regionCenter, lpIdx):
        w, c = self._getGridDimensions(randomFraction, regionSize, regionCenter, lpIdx)
        if isinstance(lensOrLensInfo, plotutil.DensInfo):
            lensInfo = lensOrLensInfo
        else:
            extraFrac = 1.0001
            lensInfo = plotutil.LensInfo(lensOrLensInfo, 
                                         bottomleft=[ c[0] - (w*extraFrac)/2, c[1] - (w*extraFrac)/2 ],
                                         topright=[ c[0] + (w*extraFrac)/2, c[1] + (w*extraFrac)/2 ])
            # Run this here already, since we know the renderer and feedbackobject here
            lensInfo.getDensityPoints(self.renderer, self.feedbackObject)
        
        return gridModule.createSubdivisionGrid(w, c, lensInfo, minSquares, maxSquares, startSubDiv)

    def setSubdivisionGrid(self, lensOrLensInfo, minSquares, maxSquares, startSubDiv = 1, randomFraction = 0.05, regionSize = None, regionCenter = None,
                           lpIdx = None):
        """Based on the lens that's provided as input, create a subdivision grid
        where regions with more mass are subdivided further, such that the number
        of resulting grid cells lies between `minSquares` and `maxSquares`. For
        more information about the procedure, see :func:`grid.createSubdivisionGrid <grale.grid.createSubdivisionGrid>`.

        The usage of the `randomFraction` parameter is the same as in :func:`setUniformGrid`.

        If specified, `regionSize` and `regionCenter` override the internally stored
        dimensions.

        In a multi-plane setting, if `lensOrLensInfo` is a :class:`MultiPlaneContainer <grale.lenses.MultiPlaneContainer>`
        result returned from a previous invert call, the specified subdivision will
        be applied to each lens plane, using the corresponding inversion result for that
        lens plane. Alternatively, a single plane lens model or :class:`LensInfo <grale.plotutil.LensInfo>`
        instance can be used for a specific lens plane, but in that case ``lpIdx`` must be
        set to the correct lens plane index.
        """
        
        if ( (not self.isMultiPlane) or
             (self.isMultiPlane and (
                isinstance(lensOrLensInfo, plotutil.DensInfo) or type(lensOrLensInfo) != lenses.MultiPlaneContainer
             ))):
            lpIdx = self._checkLensPlaneIndex(lpIdx)
            self.grid[lpIdx] = self._getSinglePlaneSubDivGrid(lensOrLensInfo, minSquares, maxSquares, startSubDiv, randomFraction, regionSize, regionCenter, lpIdx)
            return

        # Here, we have a multi-plane container, do a similar subdivision
        # step for each lens plane
        lensesAndRedshifts = lensOrLensInfo.getLensParameters()
        if len(lensesAndRedshifts) != len(self.zd):
            raise InversionException("The multi plane contain does not have the same amount of lens planes as specified for the inversion")

        for i in range(len(self.zd)):
            if not np.isclose(lensesAndRedshifts[i]["z"], self.zd[i]):
                raise InversionException("Redshift {:g} in multi plane container lens {} does not appear to match the expected redshift {:g}".format(lensesAndRedshifts[i]["z"], i+1, self.zd[i]))

            self.grid[i] = self._getSinglePlaneSubDivGrid(lensesAndRedshifts[i]["lens"], minSquares, maxSquares, startSubDiv, randomFraction, regionSize, regionCenter, i)
       
    def setDefaultInversionArguments(self, **kwargs):
        """In case you want to pass the same keyword arguments to the :func:`invert <grale.inversion.InversionWorkSpace.invert>`
        function multiple times, you can specify them using this function and omit them
        in the ``invert`` call. The different keyword arguments you can pass are the
        ones from the core :func:`invert <grale.inversion.invert>` function.
        """
        self.inversionArgs = kwargs

    def invert(self, populationSize, **kwargs):
        """For the current grid, the current images data sets, run the genetic algorithm
        for the inversion. This calls the :func:`invert <grale.inversion.invert>` function
        in this module, which you can consult for other arguments that you can specify
        using keywords. In addition to the keyword arguments specified there, you can
        also specify `gridSizeFactor`, `rescaleBasisFunctions` and `basisFunctionType`
        which will be used to set the `sizefactor`, `rescale` and `basistype` parameters
        in the `initialParameters` of the
        :func:`addBasisFunctionsBasedOnCurrentGrid <grale.inversion.InversionWorkSpace.addBasisFunctionsBasedOnCurrentGrid>`
        call that's used internally.

        If the same arguments need to be set for each call of this method, you can use
        :func:`setDefaultInversionArguments` to set them. Note that the options passed
        in ``kwargs`` will override the settings stored by that function.

        This function is now a convenience function, it first calls
        :func:`clearBasisFunctions <grale.inversion.InversionWorkSpace.clearBasisFunction>`,
        then
        :func:`addBasisFunctionsBasedOnCurrentGrid <grale.inversion.InversionWorkSpace.addBasisFunctionsBasedOnCurrentGrid>`,
        and finally
        :func:`invertBasisFunctions <grale.inversion.InversionWorkSpace.invertBasisFunctions>`.
        """
        newKwargs = { }
        newKwargs["inverter"] = self.inverter
        newKwargs["feedbackObject"] = self.feedbackObject
        for a in self.inversionArgs:
            newKwargs[a] = self.inversionArgs[a]

        for a in kwargs:
            newKwargs[a] = kwargs[a]

        self.clearBasisFunctions()

        initialParameters = { }
        if "basisFunctionType" in newKwargs:
            initialParameters["basistype"] = newKwargs["basisFunctionType"]
            del newKwargs["basisFunctionType"]

        if "rescaleBasisFunctions" in newKwargs:
            initialParameters["rescale"] = newKwargs["rescaleBasisFunctions"]
            del newKwargs["rescaleBasisFunctions"]

        if "gridSizeFactor" in newKwargs:
            initialParameters["sizefactor"] = newKwargs["gridSizeFactor"]
            del newKwargs["gridSizeFactor"]

        for i in range(len(self.zd)): # For each grid
            self.addBasisFunctionsBasedOnCurrentGrid(initialParameters=initialParameters, lpIdx=i)
        # Already processed default arguments
        return self.invertBasisFunctions(populationSize, **newKwargs, ignoreDefaultArguments=True)

    def clearBasisFunctions(self, lpIdx = "all"):
        """Clears the list of basis functions that will be used in
        :func:`invertBasisFunctions`. In a multi-plane setting
        ``lpIdx`` can be used to only clear the basis funtion for a
        specific lens plane."""
        if lpIdx == "all":
            self.basisFunctions = [ [] for i in range(len(self.zd)) ]
        else:
            lpIdx = self._checkLensPlaneIndex(lpIdx)
            self.basisFunctions[lpIdx] = [ ]

    def getBasisFunctions(self, lpIdx = None):
        """Returns the basis functions that have currently been stored, and
        of which the weights will be optimized when :func:`invertBasisFunctions`
        is called. In a multi-plane setting, ``lpIdx`` can be used to obtain
        the basis functions for a particular lens plane."""
        lpIdx = self._checkLensPlaneIndex(lpIdx)
        return self.basisFunctions[lpIdx]

    def _getImageSize(self):

        xCoords = [ ]
        yCoords = [ ]
        for imgInfo in self.imgDataList:
            
            img = imgInfo["images"]
            params = imgInfo["params"]
            imgType = params["type"]
            if imgType not in [ "extendedimages", "pointimages", "pointgroupimages" ]:
                continue
                    
            if img.getNumberOfImages() <= 1:
                continue
            
            tr = img.getTopRightCorner()
            bl = img.getBottomLeftCorner()
            xCoords.append(tr[0])
            xCoords.append(bl[0])
            yCoords.append(tr[1])
            yCoords.append(bl[1])

        if not xCoords or not yCoords:
            raise InversionException("No images are present from which the strong lensing region can be estimated")

        xMax, xMin = max(xCoords), min(xCoords)
        yMax, yMin = max(yCoords), min(yCoords)

        return ((xMax-xMin)**2 + (yMax-yMin)**2)**0.5

    def setBasisFunctions(self, basisFunctions, lpIdx=None):
        """Convenience method, just calls :func:`clearBasisFunctions` followed
        by :func:`addBasisFunctions`. In a multi-plane setting, ``lpIdx``
        needs to be used to select a particular lens plane."""
        lpIdx = self._checkLensPlaneIndex(lpIdx)
        self.clearBasisFunctions(lpIdx)
        self.addBasisFunctions(basisFunctions, lpIdx)

    def addBasisFunctions(self, basisFunctions, lpIdx=None):
        """Add basis functions that will be used in :func:`invertBasisFunctions`.
        Each entry in this list should be a dictionary with at least entries
        ``lens`` and ``center`` containing the lens model and position at which
        it should be placed. A ``mass`` parameter may be present as well,
        describing the relevant lensing mass of this basis function. If no such
        mass is specified, it will be estimated numerically as the mass inside 
        (roughly) the image region.
        
        In the final call to the :func:`invert <grale.inversion.invert>` function,
        the complete list is then passed as the `gridInfoOrBasisFunctions`
        argument. You may look there for some additional information.

        In a multi-plane setting, ``lpIdx`` needs to be used to select a particular 
        lens plane.
        """
        lpIdx = self._checkLensPlaneIndex(lpIdx)

        imgSize = None

        for bf in basisFunctions:
            d = { }
            d["lens"] = copy.deepcopy(bf["lens"])
            d["center"] = copy.deepcopy(bf["center"])
            if "mass" in bf:
                d["mass"] = copy.deepcopy(bf["mass"])
            else:
                if not imgSize:
                    imgSize = self._getImageSize() / 2.0

                lens = d["lens"]
                extra = 1.001
                lensInfo = { 
                    "lens": lens,
                    "bottomleft": [ -imgSize*extra, -imgSize*extra ],
                    "topright": [ imgSize*extra, imgSize*extra ],
                }

                mass = None
                try:
                    mass = lens.getRadialMassProfile([imgSize])[0]
                except Exception as e:
                    # TODO: for debugging
                    print(f"Unable to get mass directly: {e}")

                if mass is None:
                    # TODO: for debugging
                    print("Using numerical calculation")
                    # Abuse this function to perform the calculations
                    profile = plotutil.plotIntegratedMassProfile(lensInfo, imgSize, axes=False, 
                                                                 renderer=self.renderer,
                                                                 feedbackObject=self.feedbackObject)
                    mass = profile[1][-1] # 1 for the mass values (0 is radii), -1 for the last integrated mass
                
                d["mass"] = mass
                # TODO: for debugging
                print("Estimated mass for basis function is: {:g} solar masses".format(mass/CT.MASS_SUN))

            self.basisFunctions[lpIdx].append(d)

    def addBasisFunctionsBasedOnCurrentGrid(self, lensModelFunction = defaultLensModelFunction,
                                            initialParameters = None, lpIdx = None):
        """The goal of this function is to add basis functions based on the
        grid that's currently stored. The conversion of grid cells to basis
        functions is done using the function specified in
        `lensModelFunction`.
        
        This function takes three arguments:

         - `operation`: the name of the operation
         - `operationInfo`: information about the operation
         - `parameters`: additional parameters

        First `lensModelFunction` is called with `operation` set to ``"start"``, `operationInfo`
        set to a dictionary with entries ``workspace`` (this workspace) and ``grid`` (the current
        grid), and `parameters` set to `initialParameters`. The return value of the function, let's
        call this `lensModelFunctionParameters`, will be used as the `parameters` for subsequent
        calls to `lensModelFunction`
        
        Then, for each grid cell in the current grid, `lensModelFunction` is called with
        operation set to ``"add"``, `operationInfo` a dictionary with keys ``size`` and ``center``
        containing the length of a cell square side and (x,y) coordinates of the cell center
        respectively. The previously returned `lensModelFunctionParameters` is passed as the
        `parameters` argument.

        In this second phase, the `lensModelFunction` should return a tuple of two things:
        the lens model for that grid cell, and (optionally) the relevant lensing mass for
        this basis function. The lens model, the center and optionally this mass is then
        added to a list, which is eventually processed by the :func:`addBasisFunctions`
        procedure.

        In a multi-plane setting, ``lpIdx`` needs to be used to select a particular 
        lens plane.
        """
        lpIdx = self._checkLensPlaneIndex(lpIdx)

        tmpBasisFunctions = [ ]
        lensModelFunctionParameters = lensModelFunction("start", { 
            "grid": self.grid[lpIdx], 
            "lensplaneindex": lpIdx,
            "workspace": self }, initialParameters)

        grid = gridModule._fractionalGridToRealGrid(self.grid[lpIdx]) if type(self.grid[lpIdx]) == dict else copy.deepcopy(self.grid[lpIdx])
        for cell in grid:
            center, size = cell["center"], cell["size"]
            lens, mass = lensModelFunction("add", { "center": center, "size": size }, lensModelFunctionParameters)
            entry = { "lens": lens, "center": center }
            if mass is not None:
                entry["mass"] = mass

            tmpBasisFunctions.append(entry)

        self.addBasisFunctions(tmpBasisFunctions, lpIdx)

    def estimateStrongLensingMass(self, skipParamCheck = False):
        """Calls :func:`estimateStrongLensingMass <grale.inversion.estimateStrongLensingMass>`
        with the images and lens distance that have been set for this inversion work space."""

        # This just uses the Dd for the closest lens plane, for backward
        # compatibility
        # TODO: clean this, make this use z's only
        return estimateStrongLensingMass(self.Dd[0], self.imgDataList, skipParamCheck)

    def invertBasisFunctions(self, populationSize, **kwargs):
        """For the lensing scenario specified in this :class:`InversionWorkSpace` instance,
        this method calls the :func:`invert <grale.inversion.invert>` function, but instead
        of using a grid, a list of basis functions (arbitrary basic lens models) is used
        of which the weights need to be optimized by the genetic algorithm. 
        function.
        """
        newKwargs = { }
        newKwargs["inverter"] = self.inverter
        newKwargs["feedbackObject"] = self.feedbackObject
        newKwargs["cosmology"] = self.cosm
        
        defArgKey = "ignoreDefaultArguments"
        if not (defArgKey in kwargs and kwargs[defArgKey] == True):
            for a in self.inversionArgs:
                newKwargs[a] = self.inversionArgs[a]

        for a in kwargs:
            if a != defArgKey:
                newKwargs[a] = kwargs[a]

        if not self.isMultiPlane:
            lens = invert(self.imgDataList, self.basisFunctions[0], self.zd[0], self.Dd[0], populationSize, **newKwargs)
        else:
            bfAndZs = [ { "lenses": x[0], "z": x[1] } for x in zip(self.basisFunctions, self.zd) ]
            lens = invertMultiPlane(self.cosm, self.imgDataList, bfAndZs, populationSize, **newKwargs)
        return lens

    def calculateFitness(self, lensOrBackProjectedImages):
        """For the current grid, the current images data sets, calculate the fitness values
        for the specified parameter, which can be either a lens, or the already backprojected
        input images.
        
        When you've received the final result after an inversion,
        calling this function with that lens as its argument should yield the same fitness 
        values to a good approximation. Some differences may exist though, as different code 
        is used in the genetic algorithm.
        """
        fitnessObjectParameters = None
        moduleName = "general"
        if "fitnessObjectParameters" in self.inversionArgs:
            fitnessObjectParameters = self.inversionArgs["fitnessObjectParameters"]
        if "moduleName" in self.inversionArgs:
            moduleName = self.inversionArgs["moduleName"]

        # For a multi-plane lens, we'll back-project the images and use the
        # fitness calculation based on this.
        # Changing the image list is not really required, but it helps test
        # some multi-plane code paths
        if type(lensOrBackProjectedImages) == lenses.MultiPlaneContainer:
            newImgList = [ ]
            for i in self.imgDataList:
                params = dict(i["params"])
                params["z"] = i["z"] # TODO: this is a bit messy, clean this

                entry = dict(i) # copy
                entry["Dds"] = 0 # Set this to zero to indicate multi-plane scenario
                entry["params"] = params # use changed parameters
                newImgList.append(entry)

            lensOrBackProjectedImages = self.backProject(lensOrBackProjectedImages, None)
        else:
            newImgList = self.imgDataList

        return calculateFitness(newImgList, self.zd[0], fitnessObjectParameters, lensOrBackProjectedImages, moduleName, self.cosm)

    def backProject(self, lens, typeFilter = [ "pointimages", "extendedimages", "pointgroupimages" ]):
        """Takes the information of the images that were added using :func:`addImageDataToList`,
        and projects the points back onto the respective source planes using the lens
        model `lens`. Note that only the basic, single core procedure is used to 
        back-project the images.

        The `typeFilter` can be used to prevent that e.g. the null space grid is
        back-projected as well. If set to `None`, then all added ``ImagesData`` entries
        will be used. If it is a list of strings, then an entry will only be considered
        if the the `type` that was specified is present in this list. The `typeFilter`
        can also be a function that's used to filter the added entries: return ``True``
        to include an entry, or ``False`` to ignore it. The function will receive two
        parameters: an index into the list of ``ImagesData`` instances, and a dictionary
        that contains the ``ImagesData`` object as well as distances and parameters.
        """

        def listBasedFilter(idx, imgInfo, l):
            params = imgInfo["params"]
            if "type" in params and params["type"] in l:
                return True
            return False

        # None means no filter
        if typeFilter is None:
            filterFunction = lambda idx, imgInfo: True
        # If it's a list, it must match the 'type'
        elif type(typeFilter) == list:
            filterFunction = lambda idx, imgInfo: listBasedFilter(idx, imgInfo, typeFilter)
        # Can also be a function that takes the image position and image Info
        else:
            filterFunction = typeFilter

        if type(lens) == lenses.MultiPlaneContainer:
            # We don't actually need the grid params
            lensPlane = multiplane.MultiLensPlane(lens, [-0.001,-0.001], [0.001, 0.001], 4, 4, 
                                                  self.renderer, self.feedbackObject, self.cosm)

            def traceFunction(imgInfo, thetas):
                imgPlane = multiplane.MultiImagePlane(lensPlane, imgInfo["z"])
                return imgPlane.traceTheta(thetas)
        else:
            def traceFunction(imgInfo, thetas):
                return lens.traceTheta(imgInfo["Ds"], imgInfo["Dds"], thetas)

        bpImages = [ ]
        for i in range(len(self.imgDataList)):
            imgInfo = self.imgDataList[i]
            if filterFunction(i, copy.deepcopy(imgInfo)): # Make a copy so it can't be modified accidentally
                # Create a copy, we'll use this to store the modified points in
                # This way, the triangulation info is preserved for example.
                img = copy.deepcopy(imgInfo["images"])

                # Gather all coordinates so that the traceTheta call will be
                # somewhat more efficient
                allPoints = np.array([ img.getImagePointPosition(i, j) for i in range(img.getNumberOfImages()) for j in range(img.getNumberOfImagePoints(i)) ], dtype=np.double)
                allPoints = traceFunction(imgInfo, allPoints)

                # Save the traced points in the 'img' instance again

                offset = 0
                for i in range(img.getNumberOfImages()):
                    numImagePoints = img.getNumberOfImagePoints(i)
                    for j in range(numImagePoints):
                        img.setImagePointPosition(i, j, allPoints[offset + j])
                    offset += numImagePoints

                bpImages.append(img)

        return bpImages

    def setStrongAndWeakBasisFunctions(self, strongSubDivInfo, weakSubDiv = None, weakRegionSize = 0,
                                       weakRegionCenter = None, weakRandomFraction = "onesquare",
                                       weakMassScale = None, ignoreWLMassInMassScaleSearch = False):
        """This is a convenience function for strong&weak inversions. It sets basis functions
        based on a strong lensing grid (uniform or subdivision grid) as well as a uniform grid
        for the weak lensing region. It returns a tuple containing the created grids for strong
        and weak lensing regions respectively. Use :func:`invertBasisFunctions` to start the
        inversion.

        Arguments:
         - `strongSubDivInfo`: this can either be a number, causing a uniform subdivision grid
           to be created for the strong lensing region (see :func:`setUniformGrid`). Alternatively,
           it can be a tuple (lens, minSubDiv, maxSubDiv) which will be used to create a
           subdivision grid through a call to :func:`setSubdivisionGrid`.

        - `weakSubDiv`: this should be a number that specifies the uniform subdivision parameter
          for the :func:`setUniformGrid` call for the weak lensing region. If this is set to zero
          or ``None``, the function will only create the strong lensing part.

        - `weakRegionSize`: this specifies the size of the grid for the weak lensing area. The
          strong lensing size is the one that was specified in the constructor.

        - `weakRegionCenter`: the default is to use the same center as the strong lensing grid,
          as specified in the constructor. A different value can be set with this parameter.

        - `weakRandomFraction`: the default will add a random offset to the weak lensing grid
          that's based on the grid cell size. A different value will be passed directly to the
          :func:`setUniformGrid` function.

        - `weakMassScale`: while the strong lensing mass scale can be estimated from the observed
          multiple image sytems, a mass scale needs to be specified for the weak lensing area.
          The basis functions that are created for the uniform weak lensing grid, will have a 
          total mass that's equal to this value.

        - `ignoreWLMassInMassScaleSearch`: internally, the genetic algorithm will rescale the
          mass distribution, in some range around the estimated strong lensing mass. When weak
          lensing mass is included, this search should either be expanded (set `massScaleSearchType`
          to ``wide`` instead of ``regular``), or the masses of the basis functions in the weak
          lensing area should not be counted in this search.
        """

        if len(self.zd) != 1:
            raise InversionException("This is meant for a single lensplane")

        # Create the SL grid
        strongGrid = None

        try:
            baseLens, minDiv, maxDiv = strongSubDivInfo
            self.setSubdivisionGrid(baseLens, minDiv, maxDiv)
            strongGrid = copy.deepcopy(self.getGrid())
        except Exception as e:
            #print("DEBUG: ", e)
            pass

        if not strongGrid:
            try:
                self.setUniformGrid(strongSubDivInfo)
                strongGrid = copy.deepcopy(self.getGrid())
            except Exception as e:
                #print("DEBUG2: ", e)
                raise InversionException("'subDivInfo' must either be an integer for a uniform subdivision, or a tuple (lens, minSubDiv, maxSubDiv)")

        def dummyLensModelFunction(operation, operationInfo, parameters):
            r = defaultLensModelFunction(operation, operationInfo, parameters)
            
            # Set the mass that's counted in determining the scale factor in the GA to (almost) zero
            # (zero is not allowed by the GA)
            if operation == "add":
                return (r[0], 0.1)
            return r

        weakGrid = None
        if weakSubDiv:
            if not weakMassScale:
                raise InversionException("A mass scale for the weak lensing region is required")
            if not weakRegionSize:
                raise InversionException("A separate weak lensing region size must be specified")

            if weakRandomFraction == "onesquare":
                weakRandomFraction = 1.0/weakSubDiv

            self.setUniformGrid(weakSubDiv, randomFraction=weakRandomFraction,
                                regionSize=weakRegionSize, regionCenter=weakRegionCenter)
            weakGrid = copy.deepcopy(self.getGrid())

            lensModelFunction = dummyLensModelFunction if ignoreWLMassInMassScaleSearch else defaultLensModelFunction
            
        # Add the basisfunction
        self.clearBasisFunctions()
        self.setGrid(strongGrid)
        self.addBasisFunctionsBasedOnCurrentGrid() # TODO: options?
        if weakGrid:
            self.setGrid(weakGrid)
            self.addBasisFunctionsBasedOnCurrentGrid(lensModelFunction, { "totalmass": weakMassScale })

        return strongGrid, weakGrid

def getDefaultInverter():
    """Convenience function in this module, just calls 
    :func:`grale.inverters.getDefaultInverter`"""
    return inverters.getDefaultInverter()

def setDefaultInverter(x):
    """Convenience function in this module, just calls 
    :func:`grale.inverters.setDefaultInverter`"""
    inverters.setDefaultInverter(x)

