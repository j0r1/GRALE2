"""This module contains functions for gravitational lens inversion. The
most straightforward method is to use the class :class:`InversionWorkSpace`,
which keeps track of the settings and provides an easier interface to
the other functions. For absolute control however, the individual functions
can still be used.
"""
from . import constants as CT
from . import inverters
from . import inversionparams
from . import privutil
from . import grid as gridModule
from . import plotutil
from . import lenses
from . import multiplane
from . import paramdesc
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
            if imgType not in [ "extendedimages", "pointimages",
                                "pointgroupimages", "bayesstronglensing" ]:
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

def getDefaultModuleParameters(moduleName = "general"):
    """For the specified module name, query the default parameters."""

    return inverters.getInversionModuleDefaultConfigurationParameters(moduleName)

def getInversionModuleUsage(moduleName = "general"):
    """Returns a usage description that's provided by the specified genetic algorithm 
    module. The usage information for the ``general`` module can be viewed here: :ref:`usage <usage-module-general>`
    """
    return inverters.getInversionModuleUsage(moduleName)

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
     - `cosmology`: the cosmological model to use.
    """
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

    return inverters.calculateFitness(moduleName, inputImages, zd, fullFitnessObjParams, lens=lens, bpImages=bpImages)

def getDefaultConvergenceParameters(eaType):
    """Helper function to return the default convergence parameters for `eaType`"""
    return inversionparams.ConvergenceParameters(None, eaType).toDict()

def _adjustEAAndConvergenceParameters(eaType, eaCount, eaParams, convParams):

    paramClass = inverters._getEAParameterClass(eaType)

    # First do the EA parameters
    if eaType == "GA" or eaType == "NSGA2":
        if not "smallmutationsize" in eaParams:
            if eaCount == 1:
                eaParams["smallmutationsize"] = -1 # Large mutation setting
            else:
                # Default to small mutation settings with mutation size 0.1
                eaParams["smallmutationsize"] = 0.1*(0.5**(eaCount-2))

    # Extend with defaults and complete dictionary
    params = paramClass(**eaParams).getSettings()
    for k in params:
        eaParams[k] = params[k]

    # Then do the convergence parameters
    if eaType == "GA" or eaType == "NSGA2":
        if not "convergencefactor" in convParams:
            if eaCount == 1:
                convParams["convergencefactor"] = 0.1
            else:
                convParams["convergencefactor"] = 0.05

    # Extend again with the defaults
    defaultParams = getDefaultConvergenceParameters(eaType)
    for k in defaultParams:
        if not k in convParams:
            convParams[k] = defaultParams[k]

def getFullEASettings(eaType = "GA", geneticAlgorithmParameters = {}, convergenceParameters = {}):
    """This function is used to convert defaults for an evolutionary algorithm (EA) stage
    to full parameters. It is used by the :func:`invert` and :func:`invertMultiPlane`
    functions. The underlying code uses lists of different EA types, parameters for the
    algorithm, and parameters used to check for convergence. For convenience, certain
    predefined names can be expanded to such lists as well.

    Arguments:
    
     - `eaType`: if this is a list, it should contain entries "GA", "DE", "JADE" or
       "RND". It can also be a string, which in turn will be converted to such a list:
     
        - "GA" will be converted to [ "GA", "GA" ], representing the two stage genetic
          algorithm approach, where the first one will typically use large mutations
          and the second one smaller ones.
        - "DE" is expanded to [ "DE" ], representing a classical `diffential evolution <https://link.springer.com/article/10.1023/A:1008202821328>`_
          optimization.
        - "JADE" is interpreted as [ "JADE" ], corresponding to the `JADE <https://ieeexplore.ieee.org/document/4424751>`_
          algorithm, a self adapting differential evolution version.
        - "GA+JADE" is expanded to [ "GA", "GA", "JADE" ], the standard two-stage
          genetic algorithm followed by a JADE step to zoom in on the optimum more closely.

     - `geneticAlgorithmParameters` (should really be called just `algorithmParameters`): if
       this is a dictionary it will be replaced by a list containing as many copies of the
       dictionary as there are entries in the `eaType` list. Then, defaults for each EA
       are filled in, for settings that are not already present.

     - `convergenceParameters`: similar to the previous argument, if this is a dictionary
       it will first be copied a number of times. Then, defaults for the algorithm will be
       added.
    """

    # First, convert the eaType to a list
    if type(eaType) == list:
        allEATypes = copy.deepcopy(eaType) # Nothing to do
    else:
        if eaType == "GA":
            # This specifies the default of two GA instances in a row
            # Different mutation settings and convergence settings will
            # be set later if needed
            allEATypes = [ "GA", "GA" ]
        elif eaType == "GA+JADE":
            allEATypes = [ "GA", "GA", "JADE" ]
        elif eaType == "NSGA2":
            allEATypes = [ "NSGA2", "NSGA2" ]
        else:
            allEATypes = [ eaType ]

    numEAs = len(allEATypes)
    if type(geneticAlgorithmParameters) == list:
        allGeneticAlgorithmParameters = copy.deepcopy(geneticAlgorithmParameters) # Nothing to do
    else:
        # Make a copy of the specified parameters according to the number of EAs,
        # so that the same settings will be used as much as possible
        allGeneticAlgorithmParameters = [ copy.deepcopy(geneticAlgorithmParameters) for i in range(numEAs)]

    if type(convergenceParameters) == list:
        allConvergenceParameters = copy.deepcopy(convergenceParameters) # Nothing to do
    else:
        # Same as above
        allConvergenceParameters = [ copy.deepcopy(convergenceParameters) for i in range(numEAs) ]

    eaCounts = { }

    for ea,eaParams,convParams in zip(allEATypes, allGeneticAlgorithmParameters, allConvergenceParameters):
        if not ea in eaCounts:
            eaCounts[ea] = 0

        eaCounts[ea] += 1
        _adjustEAAndConvergenceParameters(ea, eaCounts[ea], eaParams, convParams)

    #print("allEATypes", allEATypes)
    #print("allGeneticAlgorithmParameters", allGeneticAlgorithmParameters)
    #print("allConvergenceParameters", allConvergenceParameters)
    return allEATypes, allGeneticAlgorithmParameters, allConvergenceParameters

def _invertCommon(inverter, feedbackObject, moduleName, calcType, fitnessObjectParameters,
                  massScale, DdAndZd, inputImages, getParamsFunction,
                  popSize, geneticAlgorithmParameters, returnNds, cosmology,
                  convergenceParameters, maximumGenerations,
                  multiPopulationParameters, eaType):

    allEATypes, allGeneticAlgorithmParameters, allConvergenceParameters = getFullEASettings(eaType, geneticAlgorithmParameters, convergenceParameters)
    if type(allEATypes) != list:
        raise InversionException("Expecting allEATypes to be a list")
    if type(allGeneticAlgorithmParameters) != list:
        raise InversionException("Expecting allGeneticAlgorithmParameters to be a list")
    if type(allConvergenceParameters) != list:
        raise InversionException("Expecting allConvergenceParameters to be a list")
    
    if len(allEATypes) != len(allGeneticAlgorithmParameters):
        raise InversionException("Expecting the same amount of EA types and EA parameters")
    if len(allEATypes) != len(allConvergenceParameters):
        raise InversionException("Expecting same amount of EA types and convergence parameters")

    # Set to some bad values as they don't make sense for a multi-plane inversion
    Dd, zd = DdAndZd if DdAndZd else (None, float("NaN"))

    parametricInversion = True if calcType.startswith("parametric") else False
    inverter, feedbackObject = privutil.initInverterAndFeedback(inverter, feedbackObject, parametricInversion)

    # Merge fitnessObjectParameters with defaults
    fullFitnessObjParams = _mergeModuleParameters(fitnessObjectParameters, moduleName, cosmology)

    # Force same setting of maximumGenerations if specified
    for fullConvParams in allConvergenceParameters:
        if maximumGenerations is not None:
            fullConvParams["maximumgenerations"] = maximumGenerations
    # Also use this for the MCMC settings
    for gaParams in allGeneticAlgorithmParameters:
        if maximumGenerations is not None and "samplegenerations" in gaParams:
            gaParams["samplegenerations"] = maximumGenerations

    allConvergenceParameters = [ inversionparams.ConvergenceParameters(d) for d in allConvergenceParameters ]

    multiPopParams = None if not multiPopulationParameters else inversionparams.MultiPopulationParameters(multiPopulationParameters)

    # Get massscale
    if massScale == "auto":
        massScale = estimateStrongLensingMass(Dd, inputImages, False)
    elif massScale == "auto_nocheck":
        massScale = estimateStrongLensingMass(Dd, inputImages, True)
    else:
        # Assume the mass scale is a value that needs to be used
        pass

    if massScale is None:
        feedbackObject.onStatus("Not using a mass scale")    
    else:
        feedbackObject.onStatus("Mass scale is: {:g} solar masses".format(massScale/CT.MASS_SUN))

    params = getParamsFunction(fullFitnessObjParams, massScale)

    # By setting the last parameter to None, no calculation is performed and only the
    # fitness components are returned. By doing this up front, we can log the 
    # fitness component description
    dummyFitness, fitnessComponentDescription = inverters.calculateFitness(moduleName, inputImages, zd, fullFitnessObjParams, None) 
    feedbackObject.onStatus("Fitness component order: {}".format(fitnessComponentDescription))

    result = inverter.invert(moduleName, calcType, popSize, allGeneticAlgorithmParameters, params, returnNds, allConvergenceParameters, multiPopParams, allEATypes)

    if not returnNds:
        result = (result[0], result[1], fitnessComponentDescription)
    else:
        result = (result[0], fitnessComponentDescription)

    return result

def invertMultiPlane(inputImages, basisLensesAndRedshifts, popSize, moduleName="general",
                     massScale="auto", allowNegativeValues=False,
                     baseLenses = [], sheetSearch="nosheet",
                     fitnessObjectParameters=None, massScaleSearchType="regular",
                     convergenceParameters={ }, geneticAlgorithmParameters={ },
                     returnNds=False, deviceIndex = "rotate",
                     inverter="default", feedbackObject="default",
                     maximumGenerations = None, cosmology = None,
                     multiPopulationParameters = None, eaType = "GA"):
    """Perform a multi-plane lens inversion. This is a rather low-level function,
    it may be easier to use an instance of :class:`InversionWorkSpace` instead.
    
    Arguments:

     - `inputImages`: a list of dictionaries with the following entries:
     
        - ``images``: an :class:`ImagesData<grale.images.ImagesData>` instance that
          describes the images of a source, the null space etc.
        - ``z``: the redshift of this source.
        - ``params``: for the ``general_gpu`` module, this could contain a 
          dictionary with at least a ``type`` field, which can be e.g. 
          ``pointimages`` or ``extendednullgrid`` (see the :ref:`usage <usage-module-general>` 
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

     - `baseLenses`: for each lens plane, a base lens can be specified to which corrections
       can be applied.

     - `sheetSearch`: by default, only the basis functions for the grid cells are used.
       You can also allow a mass-sheet basis function, which may be useful as this kind
       of effect is difficult to model by a grid of basis functions. To do so, set this
       to ``"genome"``, in which a mass-sheet component will be added to each lens plane.

     - `fitnessObjectParameters`: parameters for the lens inversion module for the
       generic algorithm. For the ``"general"`` module, more information can be
       found in the :ref:`usage <usage-module-general>` documentation.

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

     - `convergenceParameters`: see :func:`getFullEASettings`.

     - `geneticAlgorithmParameters`: see :func:`getFullEASettings`.

     - `returnNds`: by default, this function will return a single gravitational lens
       model. If there are several fitness measures however, the end result is actually
       a non-dominated set of models. The inversion module for the genetic algorithm has
       some default strategy for choosing one solution from this set. In case you'd like
       to get the complete non-dominated set instead, you can set this flag to ``True``.

     - `deviceIndex`: this multi-plane inversion uses a GPU to back-project the image
       data, and by setting a specific number, a specific device can be specified. To
       allow multiple GPUs to be used automatically, you can leave this to ``"rotate"``
       and use an :mod:`inverter <grale.inverters>` with as many processes as you have
       GPUs.

     - `inverter`: specifies the inverter to be used. See the :mod:`inverters<grale.inverters>`
       module for more information.

     - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

     - `maximumGenerations`: if the genetic algorithm didn't stop by itself after
       this many generations, stop it anyway. To test an inversion script completely,
       it can be useful to temporarily stop the genetic algorithm after only a small
       number of generations so that the code doesn't take long to run (this is now
       actually merged into the `convergenceParameters`)

     - `cosmology`: The cosmological model to use, to calculate the necessary angular
       diameter distances based on the specified redshifts.

     - `multiPopulationParameters`: TODO: experimental multi-population feature

     - `eaType`: see :func:`getFullEASettings`.
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
                basisLensesAndRedshifts, inputImages, baseLenses, massScale, sheetSearch,
                fullFitnessObjParams, allowNegativeValues,
                massScaleSearchType, deviceIndex)

    return _invertCommon(inverter, feedbackObject, moduleName, "multiplanegpu", fitnessObjectParameters,
                  massScale, None, inputImages, getParamsFunction, popSize,
                  geneticAlgorithmParameters, returnNds, cosmology,
                  convergenceParameters, maximumGenerations, multiPopulationParameters, eaType)

def invert(inputImages, basisFunctions, zd, Dd, popSize, moduleName = "general", massScale = "auto",
           allowNegativeValues = False, baseLens = None, 
           sheetSearch = "nosheet", fitnessObjectParameters = None, massScaleSearchType = "regular", convergenceParameters = { },
           geneticAlgorithmParameters = { }, returnNds = False, inverter = "default", feedbackObject = "default",
           cosmology = None, maximumGenerations = None, multiPopulationParameters = None, eaType = "GA",
           useImagePositionRandomization = False,
           ):
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
          ``pointimages`` or ``extendednullgrid`` (see the `usage :ref:<usage-module-general>`
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
       found in the :ref:`usage <usage-module-general>` documentation.

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

     - `convergenceParameters`: see :func:`getFullEASettings`.

     - `geneticAlgorithmParameters`: see :func:`getFullEASettings`

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

     - `maximumGenerations`: if the genetic algorithm didn't stop by itself after
       this many generations, stop it anyway. To test an inversion script completely,
       it can be useful to temporarily stop the genetic algorithm after only a small
       number of generations so that the code doesn't take long to run (this is now
       actually merged into the `convergenceParameters`)

     - `multiPopulationParameters`: TODO: experimental multi-population feature

     - `eaType`: see :func:`getFullEASettings`.
     
     - `useImagePositionRandomization`: if set to ``True``, images that contain
       the ``"positionuncertainty"`` property will be randomized slightly according
       to this property's value. This can be used to avoid overfitting.
    """

    # Check positional uncertainties
    hasPositionUncert = False
    if useImagePositionRandomization: # Check that we have something to randomize
        for i in inputImages:
            img = i["images"]
            if "positionuncertainty" in img.getKnownPropertyNames():
                hasPositionUncert = True
                break
        else:
            raise InversionException("Input image position randomization requested, but no 'positionuncertainty' property information found")
    
    initialUncertSeed = 0 if not hasPositionUncert else random.getrandbits(64)

    def getParamsFunction(fullFitnessObjParams, massScale):
        return inversionparams.LensInversionParametersSinglePlaneCPU(inputImages, basisFunctions,
                                                         Dd, zd, massScale, allowNegativeValues, baseLens, 
                                                         sheetSearch, fullFitnessObjParams, massScaleSearchType,
                                                         useImagePositionRandomization, initialUncertSeed)

    return _invertCommon(inverter, feedbackObject, moduleName, "singleplanecpu", fitnessObjectParameters,
                  massScale, [Dd, zd], inputImages, getParamsFunction, popSize,
                  geneticAlgorithmParameters, returnNds, cosmology, convergenceParameters,
                  maximumGenerations, multiPopulationParameters, eaType)

def _processCoupledParameters(varParams):

    # First analyze everything so we can find the reference variables
    cNameMap = { }

    # First keep track of which variable need to be coupled
    for x in varParams:
        if not "cname" in x:
            continue

        cname = x["cname"]

        if not cname in cNameMap:
            cNameMap[cname] = { "varparams": [] }

        cNameMap[cname]["varparams"].append(x)

    # Then find the one to use as reference
    for cname in cNameMap:
        refCandidates = []
        for p in cNameMap[cname]["varparams"]:
            if not p["ccode"]:
                refCandidates.append(p)

        if len(refCandidates) == 0:
            raise InversionException(f"No reference (unadjusted) variable found for cname '{cname}'")

        if len(refCandidates) > 1:
            print(f"WARNING: multiple reference candidates found for cname '{cname}', will be using first")
            # TODO: check if their boundary settings are the same

        ref = refCandidates[0]
        print(f"INFO: reference for cname '{cname}' is: ", ref)
        cNameMap[cname]["refparam"] = ref

    coupledVarParams = []

    mapIndices = []
    codes = []

    for idx,x in enumerate(varParams):
        if not "cname" in x:
            l = len(coupledVarParams)
            coupledVarParams.append(x)
            mapIndices.append(l)
            codes.append("")
            continue

        cname = x["cname"]
        m = cNameMap[cname]
        if m["refparam"]["name"] == x["name"]: # This is the reference one, use it
            l = len(coupledVarParams)
            coupledVarParams.append(x)
            mapIndices.append(l)
            codes.append(x["ccode"])
            m["refindex"] = l # Record this for a second pass
        else: # Don't actually include it in the coupledVarParams
            mapIndices.append(None) # placeholder, we may not know the actual index yet
            codes.append(x["ccode"])

    assert len(mapIndices) == len(varParams)

    # Second pass, fill in the 'None' values
    for idx,x in enumerate(varParams):
        if not "cname" in x:
            continue

        if mapIndices[idx] is not None:
            continue

        cname = x["cname"]
        refIdx = cNameMap[cname]["refindex"]
        mapIndices[idx] = refIdx

    originParamsMap = list(zip(mapIndices, codes))
    numOriginParams = len(coupledVarParams)

    return coupledVarParams, originParamsMap, numOriginParams

def _checkInternalRecalcLensConsistency(internalRecalcLens, parametricLensDescription):
    from .constants import ANGLE_ARCSEC, MASS_SUN
    calcLensDesc = internalRecalcLens if type(internalRecalcLens) == dict else eval(paramdesc.createParametricDescription(internalRecalcLens))

    t1, t2 = calcLensDesc["type"], parametricLensDescription["type"]
    if t1 != t2:
        raise Exception(f'Lens type mismatch for fitness calculation: {t1} != {t2}')

    p1, p2 = calcLensDesc["params"], parametricLensDescription["params"]
    if t1 == "CompositeLens":
        assert type(p1) == list and type(p2) == list, "Parameters of CompositeLens should be lists"
        assert len(p1) == len(p2), "Different number of sublenses in CompositeLens"

        for subP1, subP2 in zip(p1, p2):
            for prop in subP1:
                if not prop in subP2:
                    raise InversionException(f"Missing property {prop} in lens description")

            subLens1, subLens2 = subP1["lens"], subP2["lens"]
            _checkInternalRecalcLensConsistency(subLens1, subLens2)
    else:
        assert type(p1) == dict and type(p2) == dict, "Parameters should be contained in dictionaries"
        for k in p1:
            if not k in p2:
                raise InversionException(f"Missing property {k} in lens description")

def getDefaultsForRetraceType(typeName):
    """TODO"""

    supportedTypes = [ "NoTrace", "SingleStepNewton", "MultiStepNewton", "ExpandedMultiStepNewton" ]
    if not typeName in supportedTypes:
        raise InversionException(f"Unsupported retrace type '{typeName}', valid names are: " + ",".join(supportedTypes))

    defaults = { "type": typeName }
    if typeName == "NoTrace":
        pass # No other parameters

    elif typeName == "SingleStepNewton":
        pass # No other parameters

    elif typeName == "MultiStepNewton":
        defaults["evaluations"] = 8

    elif typeName == "ExpandedMultiStepNewton":
        defaults["evalsperpoint"] = 8
        defaults["maxgridsteps"] = 3
        defaults["acceptthreshold"] = "auto"
        defaults["gridspacing"] = "auto"
        defaults["layout"] = "FullGrid"

    return defaults

def _getScaleFromBayesStrongLensingImages(inputImages):
    allPts = []
    for img in inputImages:
        if img["params"]["type"] != "bayesstronglensing":
            continue

        imgDat = img["images"]
        for image in imgDat.getAllImagePoints():
            for x in image:
                allPts.append(x["position"])

    pts = np.array(allPts)
    maxVec = np.max(pts, axis=0)
    minVec = np.min(pts, axis=0)
    diff = maxVec-minVec
    scale = np.sum(diff**2)**0.5
    return scale

def _getAcceptThresholdFromBayesStrongLensingImages(inputImages):
    # TODO: make this configurable?
    return _getScaleFromBayesStrongLensingImages(inputImages)/10000

def _getGridSpacingFromBayesStrongLensingImages(inputImages):
    # TODO: make this configurable?
    return _getScaleFromBayesStrongLensingImages(inputImages)/100

def invertParametric(inputImages, parametricLensDescription, zd, Dd, popSize, moduleName = "general",
           defaultInitialParameterFraction = 0.1, clampToHardLimits = False, fitnessObjectParameters = None, convergenceParameters = { },
           geneticAlgorithmParameters = { }, returnNds = False, inverter = "default", feedbackObject = "default",
           cosmology = None, maximumGenerations = None, eaType = "JADE", deviceIndex = "rotate",
           useImagePositionRandomization = False, allowUnusedPriors = False,
           retraceParams = { "type": "ExpandedMultiStepNewton" },
           retraceSourcePlaneThreshold = "auto", sourcePositionEstimate = "average",
           allowEmptyInitialValueRange = False,
           forceScales = None, 
           internalRecalcLens = None, internalCalcFitnessParams = None,
           ):
    """This is a low-level function, used by the similarly named function
    in :class:`InversionWorkSpace`.

    Arguments:
     - `inputImages`: describes the observed images and their properties,
       see :func:`invert` for a full description.

     - `parametricLensDescription`: a dictionary that describes the lens model
       and which parameters of the model can vary within which bounds. See
       the :mod:`paramdesc <grale.paramdesc>` module for more information.

     - `zd`: the redshift to the lens, used when calculating time delays
       (for other purposes the angular diameter distance will be used).

     - `Dd`: the angular diameter distance to the lens.

     - `popSize`: the size of the population in the genetic algorithm, e.g. 512.

     - `moduleName`: name of the inversion module for the evolutionary algorithm.

     - `defaultInitialParameterFraction`: this is passed onto the 
       :func:`analyzeParametricLensDescription <grale.paramdesc.analyzeParametricLensDescription>`
       function.

     - `clampToHardLimits`: also passed to the same function.

     - `fitnessObjectParameters`: parameters for the lens inversion module for the
       evolutionary algorithm. For the ``"general"`` module, more information can be
       found in the :ref:`usage <usage-module-general>` documentation.

     - `convergenceParameters`: see :func:`getFullEASettings`.

     - `geneticAlgorithmParameters`: see :func:`getFullEASettings`

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

     - `maximumGenerations`: if the genetic algorithm didn't stop by itself after
       this many generations, stop it anyway. To test an inversion script completely,
       it can be useful to temporarily stop the genetic algorithm after only a small
       number of generations so that the code doesn't take long to run (this is now
       actually merged into the `convergenceParameters`)

     - `eaType`: see :func:`getFullEASettings`.

     - `deviceIndex`: this parametric inversion uses a GPU to back-project the image
       data, and by setting a specific number, a specific device can be specified. To
       allow multiple GPUs to be used automatically, you can leave this to ``"rotate"``
       and use an :mod:`inverter <grale.inverters>` with as many processes as you have
       GPUs.

     - `useImagePositionRandomization`: if set to ``True``, images that contain
       the ``"positionuncertainty"`` property will be randomized slightly according
       to this property's value. This can be used to avoid overfitting.
     
     - `allowUnusedPriors`: by default, if prior information is set on the parameters
       but the fitness measure cannot incorporate this, an error will be generated. To
       ignore this error, the flag can be set to ``True``.

     - `retraceParams`: TODO, either a dict or "disable", for the bayesian strong
       lensing

     - `sourcePositionEstimate`: TODO, "mean" or "magweighted"

     - `retraceSourcePlaneThreshold`: TODO, either "auto" or a number

     - `allowEmptyInitialValueRange`: TODO

     - `forceScales`: TODO

     - `internalRecalcLens`: For internal use

     - `internalCalcFitnessParams`: For internal use
    """

    desc = paramdesc.analyzeParametricLensDescription(parametricLensDescription, Dd, defaultInitialParameterFraction, clampToHardLimits, allowEmptyInitialValueRange,
                                                      forceScales=forceScales)
    clProbCode = desc["neglogprobcode"]
    
    templateLens = desc["templatelens"]
    deflScale, potScale = desc["scales"]["deflectionscale"], desc["scales"]["potentialscale"]
    varParams = desc["variablefloatparams"]

    # We need the original offsets
    offsets = [ x["offset"] for x in varParams ]

    isUsingCoupledParameters = sum([ 1 if "cname" in x else 0 for x in varParams ])
    originParametersMap, numOriginParams = [], 0
    if isUsingCoupledParameters:
        # override varParams to the things that are really independent
        origVarParams = varParams
        varParams, originParametersMap, numOriginParams = _processCoupledParameters(origVarParams)
    else:
        origVarParams = varParams

    initMin = [ x["initialrange"][0] for x in varParams ]
    initMax = [ x["initialrange"][1] for x in varParams ]
    hardMin = [ x["hardlimits"][0] for x in varParams ]
    hardMax = [ x["hardlimits"][1] for x in varParams ]

    if internalRecalcLens is not None:
        if internalCalcFitnessParams is not None:
            raise InversionException("Internal error: both internalCalcFitnessParams and internalRecalcLens are set")
        _checkInternalRecalcLensConsistency(internalRecalcLens, parametricLensDescription)

        _, floatParams = internalRecalcLens.getCLParameters(deflScale, potScale)
        # Change initmin/initmax to be the exact same value? Will this trigger an
        # error in a check later on?
        for idx, x in enumerate(varParams):
            off = x["offset"]
            v = floatParams[off]
            initMin[idx] = v
            initMax[idx] = v

    # TODO: should be a cleaner way to do this
    usingMcmc = False
    if eaType == "MCMC" or eaType == "MCMC-MH":
        usingMcmc = True
    elif type(eaType) == list or type(eaType) == tuple:
        for t in eaType:
            if t == "MCMC" or eaType == "MCMC-MH":
                usingMcmc = True

    infOnBoundsViolation = False if not usingMcmc else True

    if usingMcmc and useImagePositionRandomization:
        raise InversionException("Cannot use input position randomization together with MCMC, there the uncertainties need to be used in the logprob")

    # Check positional uncertainties
    hasPositionUncert = False
    if useImagePositionRandomization: # Check that we have something to randomize
        for i in inputImages:
            img = i["images"]
            if "positionuncertainty" in img.getKnownPropertyNames():
                hasPositionUncert = True
                break
        else:
            raise InversionException("Input image position randomization requested, but no 'positionuncertainty' property information found")
    
    initialUncertSeed = 0 if not hasPositionUncert else random.getrandbits(64)

    retraceImages = []
    anyRetrace = False
    for img in inputImages:
        if not "params" in img:
            raise InversionException("Input image is not a dictionary with 'params' entry")
        params = img["params"]
        if not "type" in params:
            raise InversionException("Input image parameters does not have a 'type' entry")

        # TODO: is there a better way to decide when to retrace something? For now this
        #       is only needed in this case, and can only be done by the opencl code
        retraceImages.append(True if params["type"] == "bayesstronglensing" else False)
        if retraceImages[-1]:
            anyRetrace = True

    #print("DEBUG: anyRetrace =", anyRetrace)
    #print("DEBUG: retraceImages =", retraceImages)

    if retraceParams == "disable":
        retraceParams = { "type": "NoTrace" }
        if anyRetrace:
            print(f"WARNING: disabling OpenCL based retrace code")
        for i in range(len(retraceImages)):
            retraceImages[i] = False

    if retraceSourcePlaneThreshold == "auto":
        if not anyRetrace:
            retraceSourcePlaneThreshold = -1 # Shouldn't matter, will trigger an error if something is wrong in the code
        else:
            retraceSourcePlaneThreshold = _getScaleFromBayesStrongLensingImages(inputImages)/10000 # TODO: make this configurable?

    if anyRetrace:
        print(f"DEBUG: using retraceSourcePlaneThreshold {retraceSourcePlaneThreshold/CT.ANGLE_ARCSEC} arcsec")

        mergedRetraceParams = getDefaultsForRetraceType(retraceParams["type"])
        for x in retraceParams:
            mergedRetraceParams[x] = retraceParams[x]

        if mergedRetraceParams["type"] == "ExpandedMultiStepNewton":
            if mergedRetraceParams["acceptthreshold"] == "auto":
                mergedRetraceParams["acceptthreshold"] = _getAcceptThresholdFromBayesStrongLensingImages(inputImages)
            if mergedRetraceParams["gridspacing"] == "auto":
                mergedRetraceParams["gridspacing"] = _getGridSpacingFromBayesStrongLensingImages(inputImages)

    else:
        mergedRetraceParams = getDefaultsForRetraceType("NoTrace")

    def getParamsFunction(fullFitnessObjParams, massScale):
        assert massScale is None, f"Internal error: expecting massScale to be None, but is {massScale}"
        return inversionparams.LensInversionParametersParametricSinglePlane(inputImages, Dd, zd,
                  templateLens, deflScale, potScale, offsets, initMin, initMax, hardMin, hardMax,
                  infOnBoundsViolation,
                  fullFitnessObjParams, deviceIndex,
                  useImagePositionRandomization, initialUncertSeed,
                  originParametersMap, numOriginParams, allowUnusedPriors,
                  retraceImages, mergedRetraceParams, retraceSourcePlaneThreshold,
                  clProbCode, allowEmptyInitialValueRange, internalCalcFitnessParams,
                  sourcePositionEstimate)

    results = _invertCommon(inverter, feedbackObject, moduleName, "parametricsingleplane", fitnessObjectParameters,
                  None, [Dd, zd], inputImages, getParamsFunction, popSize,
                  geneticAlgorithmParameters, returnNds, cosmology, convergenceParameters,
                  maximumGenerations, None, eaType)

    # Return both varParams and origVarParams here, the first is needed to interpret
    # the bayesian samples, the second contains information about all parameters that
    # are changed, can be useful to refine a lens description for example
    # so that info about cname is available too
    if len(results) == 2: # NDS returned
        allLenses, allFitnesses = [], []
        for x in results[0]:
            allLenses.append(x[0])
            allFitnesses.append(x[1])
        obj = { "lenses": allLenses, "fitnesses": allFitnesses,
                "fitnessorder": results[1] }
    elif len(results) == 3: # Single lens is returned
        obj = { "lens": results[0], "fitness": results[1],
                "fitnessorder": results[2] }
    else:
        raise InversionException("Internal error: unexpected number of results")

    obj["optparams"] = varParams
    obj["allvarparams"] = origVarParams
    obj["scales"] = desc["scales"]

    return obj

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
     - ``rescale`` (default is ``False``): by default, the weight of a basis function is a
       measure of its total mass. This implies that smaller grid cells will correspond
       to larger densities, i.e. if you set all weights of the basis functions to the
       same value, the regions with smaller grid cells will have a considerably larger
       density. Since the usual approach will be to subdivide the regions that contain
       more mass into smaller grid cells, this does make sense and appears to produce
       very good results in most cases. To make this effect less
       pronounced, you can set this parameter to ``True``, and the basis functions will
       be rescaled. The smaller grid cells will still have higher densities, but not
       as much as before (in case square basis functions were used, equal weights would
       generate equal densities irrespective of the grid size).
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

# We need a new type to disambiguate between a list of regions in a single
# lensplane, and regions for multiple lens planes
class Regions(object):
    """For free-form lens inversions, typically a grid based approach
    is used, where a square region is subdivided into smaller grid cells.
    A single lens plane can be covered by multiple such regions, in case
    the lens is relatively large and has distinct strong lensing parts.
    This class represents one or more of these regions, in a single
    lens plane. For a multi-plane scenario, multiple of these Regions
    objects could be used.
    """

    def __init__(self, regionInfoList):
        """Construct an instance based on `regionInfoList`. In principle
        this should be a list of dictionaries with keys ``"size"`` and
        ``"center"`` (the latter is assumed to be (0, 0) if not present).
        It can also just be a single dictionary, as a shorthand for a list
        containing only this one entry.
        """
        if type(regionInfoList) == dict:
            regionInfoList = [ regionInfoList ]

        self.regionInfoList = []
        for regionInfo in regionInfoList:
            size = regionInfo["size"]
            center = regionInfo["center"] if "center" in regionInfo else [0.0, 0.0]
            self.regionInfoList.append({ "size": size, "center": center })

    def getNumberOfRegions(self):
        """Returns the number of different regions that were specified in
        the constructor."""
        return len(self.regionInfoList)
    
    def getRegions(self):
        """Returns the regions that were specified in the constructor."""
        return self.regionInfoList

class _MultiGridWrapper(object):
    def __init__(self, singleOrMultiGrid):
        singleOrMultiGrid = copy.copy(singleOrMultiGrid)

        if type(singleOrMultiGrid) == dict:    
            singleOrMultiGrid = [ singleOrMultiGrid ]
        elif type(singleOrMultiGrid) != list:
            raise InversionException("Unexpected argument type {}".type(singleOrMultiGrid))

        newList = []
        for entry in singleOrMultiGrid:
            if not ("size" in entry and "center" in entry):
                raise Exception("No size or center present in grid/cell info")
            
            if not "cells" in entry:
                # We may be dealing with a list of absolute cells, convert to fractional
                entry["cells"] = [ { "size": 1, "center": [0,0]} ]

            newList.append(entry)

        self.multiGrid = newList            

    def getNumberOfGrids(self):
        return len(self.multiGrid)

    def getGrid(self, idx = 0):
        return copy.copy(self.multiGrid[idx])
    
    def getAllGrids(self):
        return copy.copy(self.multiGrid)

def _getLensPlaneRegions(regionSize, regionCenter, multiRegionInfo, numPlanes) -> list[gridModule.MultiGridCreator]:

    assert numPlanes >= 1, "Unexpected: numPlanes is {}".format(numPlanes)

    if numPlanes == 1:
        if multiRegionInfo is None:
            return [ gridModule.MultiGridCreator(regionSize, regionCenter) ]
        
        if type(multiRegionInfo) == Regions:
            return [ gridModule.MultiGridCreator(multiRegionInfo=multiRegionInfo.getRegions()) ]

        if type(multiRegionInfo) == dict:
            return [ gridModule.MultiGridCreator(multiRegionInfo=multiRegionInfo) ]
        
        newInfo = []
        for x in multiRegionInfo:
            if type(x) == Regions:
                newInfo.append(x.getRegions())
            elif type(x) == dict:
                newInfo.append(x)
            else:
                raise Exception("Unexpected type {} in multiRegionInfo list".format(type(x)))
        return [ gridModule.MultiGridCreator(multiRegionInfo=newInfo) ]

    # numPlanes > 1
    if multiRegionInfo is None:
        return [ gridModule.MultiGridCreator(regionSize, regionCenter) for i in range(numPlanes) ]

    if numPlanes != len(multiRegionInfo):
        raise InversionException("There must be the same amount of regions as lens planes")

    gridCreators = []
    for r in multiRegionInfo:
        if type(r) != Regions:
            raise InversionException("Each entry in multiRegionInfo must be a Regions object")
        gridCreators.append(gridModule.MultiGridCreator(multiRegionInfo=r.getRegions()))

    return gridCreators

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
    def __init__(self, zLens, regionSize = None, regionCenter = None, multiRegionInfo = None,
                 inverter = "default", 
                 renderer = "default", feedbackObject = "default", cosmology = "default"):
        
        """Constructor for this class.

        Arguments:
         - `zLens`: the redshift to the gravitational lens for a single lens plane inversion,
           or a list of redshifts for an (experimental) multi lens plane inversion.

         - `regionSize` and `regionCenter`: the width and height of the region in which the
           inversion should take plane, as well as its center. This will be used to base the
           grid dimensions on, but by default some randomness will be added (see e.g. :func:`setUniformGrid`).
           For a more complex setup, where multiple such regions (in the same lens plane) are
           needed, these can be set to `None`, and the `multiRegionInfo` parameter below can
           be used. While region information is not needed when :func:`invertBasisFunctions`
           or :func:`invertParametric` are used, at least one region needs to be defined, even
           if it is just a 'dummy' value.
           
           These are the default values in :func:`setUniformGrid` and :func:`setSubdivisionGrid`
           but can still be overridden there.

           In case of a multi-plane inversion, you can specify different sizes for each plane
           if desired, otherwise the same settings will be used for all planes.

         - `multiRegionInfo`: for a more complex lens, where several strong lensing regions are
           present, you may want to use several regions than can be subdivided as the inversion
           progresses. For a single lens plane, you can specify a :class:`Regions <grale.inversion.Regions>`
           object, a dictionary with ``"regionSize"`` and ``"regionCenter"`` keys, or a list of
           these. For a multi-plane scenario, a list of :class:`Regions <grale.inversion.Regions>`
           objects is needed, one for each lens plane.

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

        # Get the region info, store as a list of Regions, one entry in the
        # list per lensplane
        self.lensplaneRegions = _getLensPlaneRegions(regionSize, regionCenter, multiRegionInfo, len(zLens))
        assert len(self.lensplaneRegions) == len(self.zd), "There must be an equal amount of lensplanes and MultiGridCreator objects (one for each lensplane)"

        self.grid: list[_MultiGridWrapper] = [ None for z in zLens ]

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

    def setRegionSize(self, regionSize = None, regionCenter = None, 
                      multiRegionInfo = None, lpIdx = "all"):
    
        """Set the inversion region, same as in the constructor: `regionSize` and 
        `regionCenter` are the width, height and center of the region in which the
        inversion should take plane. This will be used to base the
        grid dimensions on, but by default some randomness will be added 
        (see e.g. :func:`setUniformGrid`).

        For a more complex setup, where multiple such regions (in the same lens plane) are
        needed, these can be set to `None`, and the `multiRegionInfo` parameter can
        be used. In this case, `multiRegionInfo` should be a
        :class:`Regions <grale.inversion.Regions>` object, a dictionary with
        ``"regionSize"`` and ``"regionCenter"`` keys, or a list of these.

        In a multi-plane scenario, ``lpIdx`` can be used to specify a particular
        lens plane for which the region should be set, otherwise the same region
        settings are used for each lens plane.
        """

        if lpIdx != "all":
            lpIdx = self._checkLensPlaneIndex(lpIdx)
            self.lensplaneRegions[lpIdx].setRegionSize(regionSize, regionCenter,
                                                       multiRegionInfo.getRegions() if type(multiRegionInfo) == Regions else multiRegionInfo)
        else: # all lensplanse
            for i in range(len(self.zd)):
                ri = None if not multiRegionInfo else multiRegionInfo[i]
                ri = ri.getRegions() if ri and type(ri) == Regions else ri
                self.lensplaneRegions[i].setRegionSize(regionSize, regionCenter, ri)

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
        input set. See the :ref:`usage <usage-module-general>` documentation for other types
        and parameters.

        For data with type ``"bayesellipticities"``, for the Bayesian weak lensing fitness measure,
        it is necessary to set `zs` to ``None``. That data itself contains for each individual galaxy
        the redshift and measured ellipticities, so no global redshift is allowed to be set.
        """
        # check that imgDat exists
        num = imgDat.getNumberOfImages()
        
        if zs is None:
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
        self.grid[lpIdx] = _MultiGridWrapper(grid)

    def getGrid(self, lpIdx = None):
        """Retrieves the currently set grid, e.g. for plotting using 
        :func:`plotSubdivisionGrid <grale.plotutil.plotSubdivisionGrid>`.
        
        In a multi-plane setting, ``lpIdx`` specifies a particular lens
        plane.
        """
        lpIdx = self._checkLensPlaneIndex(lpIdx)
        g = self.grid[lpIdx]
        assert type(g) == _MultiGridWrapper, "Internal error: type _MultiGridWrapper expected"
        if g.getNumberOfGrids() == 1:
            return g.getGrid(0)
        return g.getAllGrids()

    def setUniformGrid(self, subDiv, randomFraction = 0.05, regionSize = None, regionCenter = None,
                       multiRegionInfo = None, lpIdx = "all", excludeFunction = None):
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

        If specified, `regionSize`, `regionCenter` or `multiRegionInfo` override 
        the internally stored dimensions.

        In a multi-plane scenario, ``lpIdx`` can be used to specify only one
        lens plane.

        The `excludeFunction` parameter can be used to exclude certain regions
        from the generated grid, see :func:`MultiGridCreator.getUniformGrid <grale.grid.MultiGridCreator.getUniformGrid>`
        for more information.
        """
        
        def process(i) -> _MultiGridWrapper:
            mri = multiRegionInfo.getRegions() if type(multiRegionInfo) == Regions else multiRegionInfo
            mgc = gridModule.MultiGridCreator(regionSize, regionCenter, mri) if regionSize or mri else self.lensplaneRegions[i]
            assert type(mgc) == gridModule.MultiGridCreator, "Expecting mgc to be a MultiGridCreator instance"
            return _MultiGridWrapper(mgc.getUniformGrid(subDiv, randomFraction, excludeFunction))

        if lpIdx == "all":
            for i in range(len(self.grid)):
                self.grid[i] = process(i)
        else:
            lpIdx = self._checkLensPlaneIndex(lpIdx)
            self.grid[lpIdx] = process(lpIdx)

    def setSubdivisionGrid(self, lensOrLensInfo, minSquares, maxSquares, startSubDiv = 1, randomFraction = 0.05,
                           regionSize = None, regionCenter = None, multiRegionInfo = None,
                           lensFilter = None, lensInfoFilter = None, lpIdx = None, excludeFunction=None, checkSubDivFunction=None):
        """Based on the lens that's provided as input, create a subdivision grid
        where regions with more mass are subdivided further, such that the number
        of resulting grid cells lies between `minSquares` and `maxSquares`. By
        default, the procedure starts from a single large grid cell, but by
        setting `startSubDiv` to some other value N, you can start from a uniform
        NxN grid instead. For more information about the procedure,
        see :func:`MultiGridCreator.getSubdivisionGrid <grale.grid.MultiGridCreator.getSubdivisionGrid>`.

        The usage of the `randomFraction` parameter is the same as in :func:`setUniformGrid`.
        If specified, `regionSize`, `regionCenter` or `multiRegionInfo` override 
        the internally stored dimensions.

        In a multi-plane setting, if `lensOrLensInfo` is a :class:`MultiPlaneContainer <grale.lenses.MultiPlaneContainer>`
        result returned from a previous invert call, the specified subdivision will
        be applied to each lens plane, using the corresponding inversion result for that
        lens plane. Alternatively, a single plane lens model or :class:`LensInfo <grale.plotutil.LensInfo>`
        instance can be used for a specific lens plane, but in that case ``lpIdx`` must be
        set to the correct lens plane index.

        If specified, `lensFilter` can be used to apply some kind of filtering
        on the lens model, in case `lensOrLensInfo` refers to a
        :class:`gravitational lens <grale.lenses.GravitationalLens>` instance.
        Similar filtering can be done for a :class:`LensInfo <grale.plotutil.LensInfo>`
        object by specifying the `lensInfoFilter` function. The first filter function
        needs to accept four parameters: the lens object, the bottom-left
        coordinate of a region, the top-right coordinate, and an `lpIdx` parameter
        for the lens plane under consideration. For the second filter function two
        parameters are passed, the LensInfo object and an `lpIdx` parameter.
        They should return a lens or LensInfo object respectively.

        With `excludeFunction` you can determine if a grid cell should be excluded
        from the final grid, and with `checkSubDivFunction` you can specify if some
        grid cell may be subdivided further. Both parameters are passed on to
        :class:`MultiGridCreator.getSubdivisionGrid <grale.grid.MultiGridCreator.getSubdivisionGrid>`.
        """

        def process(i, lensFilter, lensInfoFilter) -> _MultiGridWrapper:
            mri = multiRegionInfo.getRegions() if type(multiRegionInfo) == Regions else multiRegionInfo
            mgc = gridModule.MultiGridCreator(regionSize, regionCenter, mri) if regionSize or multiRegionInfo else self.lensplaneRegions[i]
            assert type(mgc) == gridModule.MultiGridCreator, "Expecting mgc to be a MultiGridCreator instance"

            return _MultiGridWrapper(mgc.getSubdivisionGrid(lensOrLensInfo, minSquares, maxSquares, startSubDiv,
                                   randomFraction, lensFilter, { "lpIdx": i },
                                   lensInfoFilter, { "lpIdx": i }, excludeFunction,
                                   checkSubDivFunction, self.renderer, self.feedbackObject))

        if ( (not self.isMultiPlane) or
             (self.isMultiPlane and (
                isinstance(lensOrLensInfo, plotutil.DensInfo) or type(lensOrLensInfo) != lenses.MultiPlaneContainer
             ))):
            lpIdx = self._checkLensPlaneIndex(lpIdx)
            self.grid[lpIdx] = process(lpIdx, lensFilter, lensInfoFilter)
            return

        # Here, we have a multi-plane container, do a similar subdivision
        # step for each lens plane
        lensesAndRedshifts = lensOrLensInfo.getLensParameters()
        if len(lensesAndRedshifts) != len(self.zd):
            raise InversionException("The multi plane contain does not have the same amount of lens planes as specified for the inversion")

        for i in range(len(self.zd)):
            if not np.isclose(lensesAndRedshifts[i]["z"], self.zd[i]):
                raise InversionException("Redshift {:g} in multi plane container lens {} does not appear to match the expected redshift {:g}".format(lensesAndRedshifts[i]["z"], i+1, self.zd[i]))

            fltrLens = None if lensFilter is None else lensFilter[i]
            fltrInf = None if lensInfoFilter is None else lensInfoFilter[i]

            # TODO: use an entry from an array for excludeFunction and checkSubDivFunction?
            self.grid[i] = process(i, fltrLens, fltrInf)

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
        call that's used internally. Note that in this function call the `lensModelFunction`
        argument is set to :func:`defaultLensModelFunction <defaultLensModelFunction>`.

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

        if sum([len(bf) for bf in self.basisFunctions]) != 0:
            raise InversionException("It seems that some basisfunctions are present already, perhaps you want 'invertBasisFunctions' (or call clearBasisFunctions)")

        try:
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

            if "massScale" in newKwargs:
                initialParameters["totalmass"] = newKwargs["massScale"]

            for i in range(len(self.zd)): # For each grid
                self.addBasisFunctionsBasedOnCurrentGrid(initialParameters=initialParameters, lpIdx=i)
            # Already processed default arguments
            return self.invertBasisFunctions(populationSize, **newKwargs, ignoreDefaultArguments=True)
        finally:
            # Clear the basis functions again
            self.clearBasisFunctions()

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
        `lensModelFunction`. The default for this argument is :func:`defaultLensModelFunction <defaultLensModelFunction>`.
        
        In general, this function takes three arguments:

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

        for grid in self.grid[lpIdx].getAllGrids():
            
            tmpBasisFunctions = [ ]
            lensModelFunctionParameters = lensModelFunction("start", { 
                "grid": grid, 
                "lensplaneindex": lpIdx,
                "workspace": self }, initialParameters)

            grid = gridModule._fractionalGridToRealGrid(grid)
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
            lens = invertMultiPlane(self.imgDataList, bfAndZs, populationSize, **newKwargs)
        return lens

    def calculateParametricFitness(self, lensOrParameters, parametricLensDescription = None, returnLens = False,
                                   forceScales = None):
        """TODO"""
        # TODO: add warning that if raw parameters are specified, the parametric description
        #       needs to be exactly the same as when the parameters were created. Otherwise
        #       other deflection/potential scales could be used, causing the raw genome values
        #       to no longer correspond to the same situation

        if isinstance(lensOrParameters, lenses.GravitationalLens):
            lens = lensOrParameters
            if parametricLensDescription is None:
                from .constants import ANGLE_ARCSEC, MASS_SUN
                def cvf(value, lensName, paramName, uniqueParamName, fullParams):
                    # Consider every parameter as a variable one, but set an empty range
                    # from which this parameter can be chosen
                    return { "initmin": value, "initmax": value }

                parametricLensDescription = eval(paramdesc.createParametricDescription(lens, convertValueFunction=cvf))

            obj = self.invertParametric(parametricLensDescription, 1, eaType="CALCULATE", maximumGenerations=0, internalRecalcLens=lens, allowEmptyInitialValueRange=True, inverter="threads:1",
                                        forceScales=forceScales)
            if returnLens:
                return obj["lens"], obj["fitness"], obj["fitnessorder"]
            return obj["fitness"], obj["fitnessorder"]

        elif isinstance(lensOrParameters, np.ndarray):
            params = lensOrParameters
            if len(params.shape) == 1:
                params = params.reshape((1,-1))
                needFinalReshape = True
            elif len(params.shape) == 2:
                needFinalReshape = False
            else:
                raise InversionException("The numpy array in 'lensOrParameters' should be either 1 or 2 dimensional")

            obj = self.invertParametric(parametricLensDescription, params.shape[0], eaType="CALCULATE", maximumGenerations=0, internalCalcFitnessParams=params, returnNds=True,
                                                                              forceScales=forceScales)

            fitnesses = obj["fitnesses"]
            if needFinalReshape:
                fitnesses = fitnesses[0]

            if returnLens:
                calcLenses = [ x for x in obj["lenses"] ]
                if needFinalReshape:
                    calcLenses = calcLenses[0]
                return calcLenses, fitnesses, obj["fitnessorder"]

            return fitnesses, obj["fitnessorder"]
        else:
            raise InversionException("The 'lensOrParameters' argument must be either a GravitationalLens instanse, or a numpy array")

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
            lensPlane = multiplane.MultiLensPlane(lens, [-0.0000001,-0.0000001], [0.0000001, 0.0000001], 4, 4, 
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

    def _getStrongAndWeakGrids(self, strongSubDivInfo, weakSubDiv, weakRegionSize,
                               weakRegionCenter, weakMultiRegionInfo,
                               weakRandomFraction,
                               strongExcludeFunction, strongCheckSubDivFunction,
                               weakExcludeFunction):

        if len(self.zd) != 1:
            raise InversionException("This is meant for a single lensplane")

        # Create the SL grid
        strongGrid = None

        if strongSubDivInfo:
            
            if type(strongSubDivInfo) == dict:
                startSubDiv = 1
                baseLens, minDiv, maxDiv = strongSubDivInfo["lens"], strongSubDivInfo["mindiv"], strongSubDivInfo["maxdiv"]
                if "startsubdiv" in strongSubDivInfo:
                    startSubDiv = strongSubDivInfo["startsubdiv"]

                self.setSubdivisionGrid(baseLens, minDiv, maxDiv, startSubDiv,
                                    excludeFunction=strongExcludeFunction,
                                    checkSubDivFunction=strongCheckSubDivFunction)
            
            elif type(strongSubDivInfo) == list or type(strongSubDivInfo) == tuple: # Can either be subdiv info (old style) or uniform info
                if issubclass(type(strongSubDivInfo[0]), lenses.GravitationalLens):
                    baseLens, minDiv, maxDiv = strongSubDivInfo
    
                    self.setSubdivisionGrid(baseLens, minDiv, maxDiv, 1,
                                    excludeFunction=strongExcludeFunction,
                                    checkSubDivFunction=strongCheckSubDivFunction)
                else: # Try uniform grid, if subdiv entries for each grid
                    self.setUniformGrid(strongSubDivInfo, excludeFunction=strongExcludeFunction)
            else: # Try uniform grid, perhaps it's just a number
                self.setUniformGrid(strongSubDivInfo, excludeFunction=strongExcludeFunction)

            strongGrid = copy.deepcopy(self.getGrid())

        weakGrid = None
        if weakSubDiv:

            if weakRandomFraction == "onesquare":
                weakRandomFraction = 1.0/weakSubDiv

            self.setUniformGrid(weakSubDiv, randomFraction=weakRandomFraction,
                                regionSize=weakRegionSize, regionCenter=weakRegionCenter,
                                multiRegionInfo=weakMultiRegionInfo,
                                excludeFunction=weakExcludeFunction)
            weakGrid = copy.deepcopy(self.getGrid())

        return strongGrid, weakGrid

    def setStrongAndWeakBasisFunctions(self, strongSubDivInfo, weakSubDiv = None, weakRegionSize = None,
                                       weakRegionCenter = None,
                                       weakMultiRegionInfo = None,
                                       weakRandomFraction = "onesquare",
                                       weakMassScale = None, ignoreWLMassInMassScaleSearch = False,
                                       strongExcludeFunction = None, strongCheckSubDivFunction = None,
                                       weakExcludeFunction = None,
                                       strongLensModelFunction = defaultLensModelFunction,
                                       strongLensModelInitialParams = {},
                                       weakLensModelFunction = defaultLensModelFunction,
                                       weakLensModelInitialParams = {}
                                       ):
        """This is a convenience function for strong&weak inversions. It sets basis functions
        based on a strong lensing grid (uniform or subdivision grid) as well as a uniform grid
        for the weak lensing region. It returns a tuple containing the created grids for strong
        and weak lensing regions respectively. Use :func:`invertBasisFunctions` to start the
        inversion.

        Arguments:
         - `strongSubDivInfo`: this can either be a number, causing a uniform subdivision grid
           to be created for the strong lensing region (see :func:`setUniformGrid`). Alternatively,
           it can be a tuple (lens, minSubDiv, maxSubDiv) which will be used to create a
           subdivision grid through a call to :func:`setSubdivisionGrid`, or more explicitly
           a dictionary containing keys ``"lens"``, ``"mindiv"``, ``"maxdiv"`` and optionally
           ``"startsubdiv"``.

        - `weakSubDiv`: this should be a number that specifies the uniform subdivision parameter
          for the :func:`setUniformGrid` call for the weak lensing region. If this is set to zero
          or ``None``, the function will only create the strong lensing part.

        - `weakRegionSize` and `weakRegionCenter`: these specify the dimensions of the grid for 
          the weak lensing area. The strong lensing size is by default the one that was specified
          in the constructor. If the weak lensing part should consist of different regions, 
          the following `weakMultiRegionInfo` parameter can be used.
        
        - `weakMultiRegionInfo`: similar to `multiRegionInfo` from the constructor, but for the
          weak lensing grid. This allows more complex regions to be specified.

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

        - `strongExcludeFunction` and `strongCheckSubDivFunction`: for the strong lensing part,
          these are passed as `excludeFunction` and `checkSubDivFunction` arguments of
          :func:`setUniformGrid` and :func:`setSubdivisionGrid`.

        - `weakExcludeFunction`: for the weak lensing part, this is passed as the `excludeFunction`
          argument of :func:`setUniformGrid`.

        - `strongLensModelFunction` and `strongLensModelInitialParams`: based on the strong lensing
          grid that's built, :func:`addBasisFunctionsBasedOnCurrentGrid` will be called to actually
          add basis functions. These parameters are passed as `lensModelFunction` and
          `lensModelFunctionParameters` of that function respectively.

        - `weakLensModelFunction` and `weakLensModelInitialParams`: same, but for the weak lensing
          grid.
        """

        strongGrid, weakGrid = self._getStrongAndWeakGrids(strongSubDivInfo, weakSubDiv, weakRegionSize,
                                                           weakRegionCenter, weakMultiRegionInfo,
                                                           weakRandomFraction,
                                                           strongExcludeFunction, strongCheckSubDivFunction,
                                                           weakExcludeFunction)

        def dummyLensModelFunction(operation, operationInfo, parameters):
            r = weakLensModelFunction(operation, operationInfo, parameters)
            
            # Set the mass that's counted in determining the scale factor in the GA to (almost) zero
            # (zero is not allowed by the GA)
            if operation == "add":
                return (r[0], 0.1)
            return r

        # Add the basisfunction
        self.clearBasisFunctions()
        if strongGrid:
            self.setGrid(strongGrid)
            self.addBasisFunctionsBasedOnCurrentGrid(strongLensModelFunction, strongLensModelInitialParams)

        if weakGrid:
            self.setGrid(weakGrid)

            if not weakMassScale:
                raise InversionException("A mass scale for the weak lensing region is required")

            lensModelFunction = dummyLensModelFunction if ignoreWLMassInMassScaleSearch else weakLensModelFunction
            params = { "totalmass": weakMassScale }
            for k in weakLensModelInitialParams:
                params[k] = weakLensModelInitialParams[k]
            self.addBasisFunctionsBasedOnCurrentGrid(lensModelFunction, params)

        if not strongGrid and not weakGrid:
            raise InversionException("Neither a grid for the strong lensing regions(s) nor for the weak lensing one was specified")

        return strongGrid, weakGrid
    
    def invertParametric(self, parametricLensDescription, populationSize, **kwargs):
        """This function uses the images that were added to the workspace to
        perform a parametric inversion using the :func:`invertParametric` function.
        A lens description and population size must be passed, other options to
        that function can be specified as keyword arguments.
        """

        zd = self.zd[0]
        Dd = self.Dd[0]
        imgList = self.imgDataList

        newKwargs = { }
        newKwargs["inverter"] = self.inverter
        newKwargs["feedbackObject"] = self.feedbackObject
        newKwargs["cosmology"] = self.cosm

        for a in self.inversionArgs:
            newKwargs[a] = self.inversionArgs[a]

        for a in kwargs:
            newKwargs[a] = kwargs[a]

        return invertParametric(imgList, parametricLensDescription, zd, Dd, populationSize, **newKwargs)

def getDefaultInverter():
    """Convenience function in this module, just calls 
    :func:`grale.inverters.getDefaultInverter`"""
    return inverters.getDefaultInverter()

def setDefaultInverter(x):
    """Convenience function in this module, just calls 
    :func:`grale.inverters.setDefaultInverter`"""
    inverters.setDefaultInverter(x)

