"""This module contains tools for parametric inversion: something to analyze
a description of a lens model that can be optimized parametrically, and
a routine to start such a description based on an existing lens model."""

from .constants import *
from . import lenses
import pprint
import copy
import math
import uuid
import json

class ParametricDescriptionException(Exception):
    """An exception that will be thrown in case something goes wrong when
    analyzing or creating a lens description for parametric inversion."""
    pass

# A fixed value van also be { "fixed": value } or { "fixed": value, "varname": "x" },
# to be able to assign a variable name to a fixed value (to be able to use this in
# OpenCL prior code
def _isFixedDictValue(d):
    if not "fixed" in d:
        return False

    dc = d.copy()
    del dc["fixed"]

    if "varname" in dc:
        del dc["varname"]

    if dc:
        raise ParametricDescriptionException(f"Unused keys {list(dc.keys())} in fixed value specification {d}")
    return True

def _getInitialParameterValue(params, paramKey):
    value = params[paramKey]
    if type(value) == dict:
        if _isFixedDictValue(value):
            return value["fixed"], value["varname"] if "varname" in value else True

        return (value["initmin"] + value["initmax"])*0.5, False
    
    if type(value) == list or type(value) == tuple:
        if len(value) == 1 or len(value) == 2 or len(value) == 3:
            return value[0], False
        raise ParametricDescriptionException("Too many entries for parameter")
    
    return value, True

def _getInitialMinOrMaxParameterValue(params, paramKey, fraction, isMin):
    value = params[paramKey]

    key = "initmin" if isMin else "initmax"
    fracSign = -1 if isMin else 1

    if type(value) == dict:
        if _isFixedDictValue(value):
            return value["fixed"], value["varname"] if "varname" in value else True

        return value[key], False
    
    if type(value) == list or type(value) == tuple:

        if len(value) < 1 or len(value) > 3:
            raise ParametricDescriptionException("Incorrect number of entries for parameter")

        if value[0] < 0:
            fracSign = -fracSign

        if len(value) == 1:
            fullFrac = 1 + fracSign*abs(fraction)
        else: # 2 params
            fullFrac = 1 + fracSign*abs(value[1])
        
        return value[0]*fullFrac, False
        
    return value, True

def _checkParameterValues(lensParams, paramMapping, getParamValue, positionNames, fixedPositionNames):
    remainingLensParams = lensParams.copy() # shallow copy is enough
    newParams = { }
    
    for x,y in paramMapping:
        
        newParams[x], isFixedOrName = getParamValue(lensParams, x)
        del remainingLensParams[x]

        if isFixedOrName is False:
            positionNames.append(y)
        else:
            fixedPositionNames.append(y)

    return remainingLensParams, newParams

def _processPlummerLens(lensParams, Dd, getParamValue, saveCLCode):
    positionNames, fixedPositionNames = [ ], [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("mass", "mass_scaled"),
                                                                ("width", "width_scaled")],
                                                  getParamValue, positionNames, fixedPositionNames)
    if lensParams:
        raise ParametricDescriptionException("Excess parameters for PlummerLens")
    
    return newParams, positionNames, fixedPositionNames

def _processCompositeLens(lensParams, Dd, getParamValue, saveCLCode):
    compParams = []
    positionNames, fixedPositionNames = [ ], [ ]

    for idx, subParams in enumerate(lensParams):
        subParams, newParams = _checkParameterValues(subParams, 
                                                     [ ("x", f"x_{idx}_scaled"),
                                                       ("y", f"y_{idx}_scaled"),
                                                       ("factor", f"factor_{idx}"),
                                                       ("angle", f"angle_{idx}") ],
                                                     getParamValue, positionNames, fixedPositionNames)

        l, posNames, fixedPosNames = _createTemplateLens_helper(subParams["lens"], Dd, getParamValue, saveCLCode)
        newParams["lens"] = l
        del subParams["lens"]
        for p in posNames:
            positionNames.append(f"lens_{idx}," + p)
        for p in fixedPosNames:
            fixedPositionNames.append(f"lens_{idx}," + p)
        
        if subParams:
            raise ParametricDescriptionException("Excess parameters for CompositeLens")
        
        compParams.append(newParams)
        
    return compParams, positionNames, fixedPositionNames

def _processNSIELens(lensParams, Dd, getParamValue, saveCLCode):
    positionNames, fixedPositionNames = [ ], [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("ellipticity", "ellipticity"),
                                                                ("coreRadius", "core_scaled"),
                                                                ("velocityDispersion", "sigma_scaled") ],
                                                  getParamValue, positionNames, fixedPositionNames)
    if lensParams:
        raise ParametricDescriptionException("Excess parameters for NSIELens")
    
    return newParams, positionNames, fixedPositionNames

def _processSISLens(lensParams, Dd, getParamValue, saveCLCode):
    positionNames, fixedPositionNames = [ ], [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("velocityDispersion", "sigma_scaled") ],
                                                  getParamValue, positionNames, fixedPositionNames)
    if lensParams:
        raise ParametricDescriptionException("Excess parameters for SISLens")
    
    return newParams, positionNames, fixedPositionNames

def _processMassSheetLens(lensParams, Dd, getParamValue, saveCLCode):
    positionNames, fixedPositionNames = [ ], [ ]
    if "Ds" in lensParams or "Dds" in lensParams:
        raise ParametricDescriptionException("For a mass sheet lens, 'Ds' and 'Dds' cannot be used, use 'density' instead")

    lensParams, newParams = _checkParameterValues(lensParams, [ ("density", "density_scaled") ],
                                                  getParamValue, positionNames, fixedPositionNames)

    if lensParams:
        raise ParametricDescriptionException("Excess parameters for MassSheetLens")
    return newParams, positionNames, fixedPositionNames

def _processMultiplePlummerLens(lensParams, Dd, getParamValue, saveCLCode):
    positionNames, fixedPositionNames = [ ], [ ]
    newMultiParams = []
    for i,subParams in enumerate(lensParams):
        subParams, newParams = _checkParameterValues(subParams, [ ("x", f"x_{i}_scaled"),
                                                                  ("y", f"y_{i}_scaled"),
                                                                  ("mass", f"mass_{i}_scaled"),
                                                                  ("width", f"width_{i}_scaled"),
                                                                  ],
                                                     getParamValue, positionNames, fixedPositionNames)
        if subParams:
            raise ParametricDescriptionException("Excess parameters for MultiplePlummerLens")
        
        newMultiParams.append(newParams)

    return newMultiParams, positionNames, fixedPositionNames

def _processDeflectionGridLens(lensParams, Dd, getParamValue, saveCLCode):
    # No parameters can be changed
    lpCopy = lensParams.copy()
    for k in [ "bottomleft", "topright" ]:
        if not k in lpCopy:
            raise ParametricDescriptionException(f"Expecting key '{k}' in DeflectionGridLens parameters")

        try:
            x, y = lpCopy[k][0], lpCopy[k][1]
            x, y = float(x), float(y)
        except Exception as e:
            raise ParametricDescriptionException(f"Value for '{k}' should be a fixed 2D coordinate") from e

        del lpCopy[k]

    if not "angles" in lpCopy:
        raise ParametricDescriptionException("Expecting key 'angles' in DeflectionGridLens parameters")

    angles = lpCopy["angles"]
    try:
        shp = angles.shape
    except Exception as e:
        raise ParametricDescriptionException("Expeting the 'angles' in DeflectionGridLens parameters to be a fixed numpy 2D array")

    del lpCopy["angles"]
    if lpCopy:
        raise ParametricDescriptionException("Excess parameters for DeflectionGridLens")

    return lensParams.copy(), [], [] # No variable parameters or fixed parameters that can be referred to

def _processPIMDLens(lensParams, Dd, getParamValue, saveCLCode):

    positionNames, fixedPositionNames = [ ], [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("coreradius", "coreradius_scaled"),
                                                                ("scaleradius", "scaleradius_scaled"),
                                                                ("centraldensity", "centraldensity_scaled")],
                                                  getParamValue, positionNames, fixedPositionNames)
    if lensParams:
        raise ParametricDescriptionException("Excess parameters for PIMDLens")

    return newParams, positionNames, fixedPositionNames

def _processPIEMDLens(lensParams, Dd, getParamValue, saveCLCode):

    positionNames, fixedPositionNames = [ ], [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("epsilon", "epsilon"),
                                                                ("coreradius", "coreradius_scaled"),
                                                                ("scaleradius", "scaleradius_scaled"),
                                                                ("centraldensity", "centraldensity_scaled")],
                                                  getParamValue, positionNames, fixedPositionNames)
    if lensParams:
        raise ParametricDescriptionException("Excess parameters for PIEMDLens")

    return newParams, positionNames, fixedPositionNames

def _processLTPIEMDLens(lensParams, Dd, getParamValue, saveCLCode):

    positionNames, fixedPositionNames = [ ], [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("ellipticity", "ellipticity"),
                                                                ("coreradius", "coreradius_scaled"),
                                                                ("scaleradius", "scaleradius_scaled"),
                                                                ("velocitydispersion", "velocitydispersion_scaled")],
                                                  getParamValue, positionNames, fixedPositionNames)
    if lensParams:
        raise ParametricDescriptionException("Excess parameters for LTPIEMDLens")

    return newParams, positionNames, fixedPositionNames

def _processLTPIMDLens(lensParams, Dd, getParamValue, saveCLCode):

    positionNames, fixedPositionNames = [ ], [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("coreradius", "coreradius_scaled"),
                                                                ("scaleradius", "scaleradius_scaled"),
                                                                ("velocitydispersion", "velocitydispersion_scaled")],
                                                  getParamValue, positionNames, fixedPositionNames)
    if lensParams:
        raise ParametricDescriptionException("Excess parameters for LTPIMDLens")

    return newParams, positionNames, fixedPositionNames

def _createTemplateLens_helper(parametricLensDescription, Dd, getParamValue, saveCLCode = None):
    lensType = parametricLensDescription["type"]
    lensParams = parametricLensDescription["params"]
    if not lensType in _supportedLensTypes:
        raise ParametricDescriptionException("Unknown lens type for parametric inversion")
    
    t = _supportedLensTypes[lensType]
    handler, lensClass = t["handler"], t["lens"]

    newParams, positionNames, fixedPositionNames = handler(lensParams, Dd, getParamValue, saveCLCode) # Need Dd as a parameter because of possible recursion
    lens = lensClass(Dd, newParams)

    if "neglogprob" in parametricLensDescription and saveCLCode is not None:
        saveCLCode(parametricLensDescription["neglogprob"])

    return lens, positionNames, fixedPositionNames

def _getCouplingNameAndCode(cname:str):
    if not ";" in cname:
        return cname.strip(), ""

    idx = cname.find(";")
    cn, code = cname[:idx].strip(), cname[idx+1:].strip()

    if not code.startswith("{"): # need to wrap it
        code = "{ return " + code + "; }"

    return cn, code

def _createParamOffsetInfo(l, paramNames, couplingNames, tagNames, varNames, priors, fixedParamNames, fixedVarNames,
                           deflectionscale, potentialscale):

    assert len(paramNames) == len(couplingNames), "Internal error: parameter name array and coupling name array should have same length"
    assert len(paramNames) == len(tagNames), "Internal error: parameter name array and tag name array differ in length"
    assert len(paramNames) == len(varNames), "Internal error: parameter name array and var name array differ in length"
    assert len(priors) == len(paramNames), "Internal error: parameter name array and prior array should have same length"

    assert len(fixedParamNames) == len(fixedVarNames), "Internal error: expecting same amount of fixed var names as parameters"

    adjustableParams = l.getCLAdjustableFloatingPointParameterInfo(deflectionscale, potentialscale)
    adjustableParamsDict = { }
    for p in adjustableParams:
        name = p["name"]
        if name in adjustableParamsDict:
            raise ParametricDescriptionException(f"Internal error: name {name} already exists in adjustable params dict")
        adjustableParamsDict[name] = p

    paramOffsetInfo = []    
    for n, cn, t, var, pr in zip(paramNames, couplingNames, tagNames, varNames, priors):
        if not n in adjustableParamsDict:
            raise ParametricDescriptionException(f"Internal error: specified adjustable param {n} does not seem to be valid")
        d = adjustableParamsDict[n].copy()
        if cn:
            d["cname_orig"] = cn
            d["cname"], d["ccode"] = _getCouplingNameAndCode(cn)
        if t:
            d["tag"] = t
        if var:
            d["varname"] = var
        if pr:
            _checkPrior(pr)
            d["prior"] = copy.deepcopy(pr)

        paramOffsetInfo.append(d)

    fixedOffsetInfo = []
    for n, var in zip(fixedParamNames, fixedVarNames):
        if not n in adjustableParamsDict:
            raise ParametricDescriptionException(f"Internal error: specified adjustable param {n} does not seem to be valid (for fixed param)")

        d = adjustableParamsDict[n].copy()
        if var:
            d["varname"] = var

        del d["hard_max"]
        del d["hard_min"]
        fixedOffsetInfo.append(d)

    return sorted(paramOffsetInfo, key=lambda x: x["offset"]), sorted(fixedOffsetInfo, key=lambda x: x["offset"])

def _createTemplateLens(parametricLensDescription, Dd):

    couplingNames = []
    tagNames = []
    priors = []
    varNames = []
    fixedVarNames = []

    # This is probably not the cleanest way
    def recordParamCouplingNamesWrapper(params, paramKey):
        value, isFixedOrName = _getInitialParameterValue(params, paramKey)
        if isFixedOrName is False:
            v = params[paramKey]
            n, t, pr, var = None, None, None, None
            if type(v) == dict:
                if "cname" in v:
                    n = v["cname"]
                if "tag" in v:
                    t = v["tag"]
                if "prior" in v:
                    pr = copy.deepcopy(v["prior"])
                if "varname" in v:
                    var = v["varname"]

            couplingNames.append(n)
            tagNames.append(t)
            priors.append(pr)
            varNames.append(var)

        else: # Fixed value, possibly with a name
            if type(isFixedOrName) != bool:
                assert type(isFixedOrName) == str, f"Internal error: expecting a 'str', but got '{type(isFixedOrName)}'"
                fixedVarNames.append(isFixedOrName)
            else:
                assert isFixedOrName is True, "Internal error: expecting a 'True' value here"
                fixedVarNames.append(None)

        return value, isFixedOrName

    allClCodeInfo = []
    def saveCLCode(p):
        if type(p) is list:
            for part in p:
                saveCLCode(part)
        elif type(p) is dict:
            if "clcall" in p and "clcode" in p:
                raise ParametricDescriptionException(f"Both 'clcode' and 'clcall' are present in 'clcode' {p}")

            cp = copy.copy(p)

            if "clcall" in p:
                del cp["clcall"]

                if not "args" in cp:
                    raise ParametricDescriptionException(f"Need an 'args' section in 'neglogprob' entry {p}")
                del cp["args"]
            else:
                if not "clcode" in p:
                    raise ParametricDescriptionException(f"The 'neglogprob' section doesn't contain a 'clcode' entry or 'clcall' entry: {p}")

                del cp["clcode"]

                # Can exist without args, to introduce a new helper function
                if "args" in cp:
                    del cp["args"]

            if cp:
                raise ParametricDescriptionException(f"The 'neglogprob' section contains more entries than 'args' and 'clcode'/'clcall': {list(cp.keys())}")

            allClCodeInfo.append(p)
        else:
            raise ParametricDescriptionException(f"Invalid 'neglogprob' parameter type for '{p}'")

    l, paramNames, fixedParamNames = _createTemplateLens_helper(parametricLensDescription, Dd, recordParamCouplingNamesWrapper, saveCLCode)

    scales = l.getSuggestedScales()
    intParam, floatParams = l.getCLParameters(**scales)
    varParamOffsetInfo, fixedParamOffsetInfo = _createParamOffsetInfo(l, paramNames, couplingNames, tagNames, varNames, priors,
                                                                      fixedParamNames, fixedVarNames, **scales)

    ret = { "templatelens": l,
            "paramoffsets": varParamOffsetInfo,
            "paramnames": paramNames, # in original order
            "fixedparamoffsets": fixedParamOffsetInfo,
            "scales": scales,
            "floatparams": floatParams,
            "description": copy.deepcopy(parametricLensDescription),
            "allclcodeinfo": allClCodeInfo
          }
    return ret

def _getShouldBeSameArray(templateLensDescription):
    numParams = templateLensDescription["floatparams"].shape[0]
    shouldBeSame = [ True for _ in range(numParams) ]
    for paramInfo in templateLensDescription["paramoffsets"]:
        offset = paramInfo["offset"]
        assert shouldBeSame[offset], "Internal error: changing same floating point value twice"

        shouldBeSame[offset] = False

    return shouldBeSame

def _checkShouldBeSame(templateLensDescription, minParams, maxParams, allowSameMinMax = False):
    shouldBeSame = _getShouldBeSameArray(templateLensDescription)

    for idx, value in enumerate(shouldBeSame):
        # print(f"Comparing index {idx}: {initMinParams[idx]} vs {initMaxParams[idx]}")
        if value: # everything should be same
            assert minParams[idx] == maxParams[idx], "Internal error: min/max should be same at this position"
            assert minParams[idx] == templateLensDescription["floatparams"][idx], "Internal error: min and template values should match at this position"

        else:
            if minParams[idx] == maxParams[idx]:
                if not allowSameMinMax:
                    raise ParametricDescriptionException(f"Values at floating point offset {idx} should be allowed to change, but initial min/max values are the same (so no variation will be introduced)")

            if minParams[idx] > maxParams[idx]:
                raise ParametricDescriptionException(f"Min/max value at floating point offset {idx} should be other way around")

def _createInitialMinMaxParameters(templateLensDesciption, defaultFraction, allowSameMinMax):
    
    scales = templateLensDesciption["scales"]
    Dd = templateLensDesciption["templatelens"].getLensDistance()
    parametricLensDescription = templateLensDesciption["description"]
    
    l, paramNames, _ = _createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: _getInitialMinOrMaxParameterValue(params, key, defaultFraction, True))
    _, initMinParams = l.getCLParameters(**scales)

    l, paramNames, _ = _createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: _getInitialMinOrMaxParameterValue(params, key, defaultFraction, False))
    _, initMaxParams = l.getCLParameters(**scales)

    assert templateLensDesciption["floatparams"].shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes"
    assert initMaxParams.shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes (2)"

    _checkShouldBeSame(templateLensDesciption, initMinParams, initMaxParams, allowSameMinMax)

    return initMinParams, initMaxParams

def _getHardMinOrMaxParameterValue(params, paramKey, isMin, knownParamNames, hardInfinite):

    # We're expecting either a fixed value, or a dict { "start", "hardmin", "hardmax"}
    value = params[paramKey]
    if type(value) != dict:
        return value, True

    if _isFixedDictValue(value):
        return value["fixed"], value["varname"] if "varname" in value else True

    paramName = knownParamNames.pop(0)

    key, infValue = ("hardmin", float("-inf")) if isMin else ("hardmax", float("inf"))
    extValue = value[key]
    if math.isinf(extValue):
        paramValue = value["start"] # Don't use inf as a lens parameter
        if extValue != infValue:
            raise ParametricDescriptionException(f"Got {extValue} but expecting {infValue} for {paramName}")
        hardInfinite.append(paramName) # Remember this parameter name
    else:
        paramValue = extValue

    return paramValue, False
    
def _mergeHardMinOrMaxParameterValue(params, key, knownParamNames, paramOffsets):

    # print("knownParamNames")
    # pprint.pprint(knownParamNames)
    # print("paramOffsets")
    # pprint.pprint(paramOffsets)
    startValue, isFixedOrName = _getInitialParameterValue(params, key)
    if isFixedOrName is not False:
        return startValue, isFixedOrName

    paramName = knownParamNames.pop(0)
    paramInfo = paramOffsets[paramName]

    value = params[key]

    newDict = { "start": startValue }

    if type(value) == dict:
        for pyKey, cKey in [ ("hardmin", "hard_min"), ("hardmax", "hard_max")]:
            if pyKey in value:
                newDict[pyKey] = value[pyKey]
            else:
                newDict[pyKey] = paramInfo[cKey]

    elif (type(value) == list or type(value) == tuple) and len(value) == 3:
        # Third is a percentage describing the hard bounds
        minVal, maxVal = startValue*(1-value[2]), startValue*(1+value[2])
        if startValue < 0:
            minVal, maxVal = maxVal, minVal
        newDict["hardmin"] = minVal
        newDict["hardmax"] = maxVal
    else:
        for pyKey, cKey in [ ("hardmin", "hard_min"), ("hardmax", "hard_max")]:
            newDict[pyKey] = paramInfo[cKey]
    
    # Change it!
    params[key] = newDict
    return startValue, False

def _createHardMinMaxParameters(templateLensDesciption):

    # We're going to merge this with hard min/max values
    parametricLensDescription = copy.deepcopy(templateLensDesciption["description"])
    Dd = templateLensDesciption["templatelens"].getLensDistance()
    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    paramOffsets = templateLensDesciption["paramoffsets"]
    paramOffsets = { x["name"]:x for x in paramOffsets }
    scales = templateLensDesciption["scales"]

    l, paramNames, _ = _createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: _mergeHardMinOrMaxParameterValue(params, key, knownParamNames, paramOffsets))

    # Now parametricLensDescription is modified, to contain
    # "start", "hardmin" and "hardmax" values for all variable
    # entries

    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    hardInfMinNames = []
    l, _, _ = _createTemplateLens_helper(parametricLensDescription, Dd,
                                     lambda params, key: _getHardMinOrMaxParameterValue(params, key, True, knownParamNames, hardInfMinNames))

    _, hardMinParams = l.getCLParameters(**scales) # These should all be finite


    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    hardInfMaxNames = []
    l, _, _ = _createTemplateLens_helper(parametricLensDescription, Dd,
                                     lambda params, key: _getHardMinOrMaxParameterValue(params, key, False, knownParamNames, hardInfMaxNames))

    _, hardMaxParams = l.getCLParameters(**scales) # These should all be finite

    def isNumber(num):
        return not (math.isinf(num) or math.isnan(num))

    # print(templateLensDesciption["paramnames"])
    # pprint.pprint(templateLensDesciption["floatparams"])
    # pprint.pprint(hardMinParams)
    # pprint.pprint(hardMaxParams)

    for i in range(hardMinParams.shape[0]):
        assert isNumber(hardMinParams[i]), f"Internal error: expecting all hard min parameters to start as a real number ({hardMinParams[i]})"

    for i in range(hardMaxParams.shape[0]):
        assert isNumber(hardMaxParams[i]), f"Internal error: expecting all hard max parameters to start as a real number ({hardMaxParams[i]})"
    
    for paramName in hardInfMinNames:
        offset = paramOffsets[paramName]["offset"]
        hardMinParams[offset] = float("-inf")

    for paramName in hardInfMaxNames:
        offset = paramOffsets[paramName]["offset"]
        hardMaxParams[offset] = float("inf")

    assert templateLensDesciption["floatparams"].shape == hardMinParams.shape, "Internal error: mismatch in floating point parameter shapes (3)"
    assert hardMaxParams.shape == hardMinParams.shape, "Internal error: mismatch in floating point parameter shapes (4)"

    _checkShouldBeSame(templateLensDesciption, hardMinParams, hardMaxParams)

    return hardMinParams, hardMaxParams

def analyzeParametricLensDescription(parametricLens, Dd, defaultFraction, clampToHardLimits = False, allowSameInitialMinMax = False):
    """Analyze the parametric lens description in `parametricLens`, which
    should be a dictionary, for a lens at angular diameter distance `Dd`.
    In case a parameter is set to change with some fraction about a value,
    `defaultFraction` is used if no specific fraction is specified. By
    default, an error will be generated if some initial value bounds exceed the
    hard bounds, but if `clampToHardLimits` is set to ``True``, the hard
    limit will be used as bound for the initial value.
    
    Returns a dictionary with the following entries:

     - ``templatelens``: a lens model constructed from the description
     - ``floatparams``: the floating point parameters for the model, some of
       which can be changed.
     - ``scales``: itself a dictionary with entries for the angular and potential
       scales that are used.
     - ``variablefloatparams``: describes which of the floating point parameters
       can be changed, and within which boundaries.

    An example `parametricLens` object:
    
    .. code:: python

       parametricLens = {
         "type": "MultiplePlummerLens",
         "params": [
            { 
               # This represents a fixed value 
               "x": 2*ANGLE_ARCSEC, 
               # This specifies an initial value that can change within 20% of the specified one.
               # During optimization these bounds can be exceeded
               "y": [ 1*ANGLE_ARCSEC, 0.20 ] 
               # Also a variable that can change, the amount determined by `defaultFraction`
               "mass": [ 1e13*MASS_SUN ],
               # When three entries are specified, the second is again to determine the
               # range of initial values, the last one is a fraction to fix hard limits.
               "width": [ 3*ANGLE_ARCSEC, 0.1, 0.5 ]
            },
            {
                # Here the initial values will be chosen in the specified interval
                "x": { "initmin": -2*ANGLE_ARCSEC, "initmax": 2*ANGLE_ARCSEC },
                # Similar, but hard bounds for the parameter are specified as well
                # In case none are specified, default values are used for these,
                # depending on the parameter type. Specifying them overrides the
                # defaults.
                "y": { "initmin": -2*ANGLE_ARCSEC, "initmax": 2*ANGLE_ARCSEC,
                       "hardmin": -3*ANGLE_ARCSEC, "hardmax": 3*ANGLE_ARCSEC },
                
                # And some fixed values
                "mass": 1e13*MASS_SUN,
                "width": 2*ANGLE_ARCSEC
            }
         ]
       }

    When a dictionary is used to specify the initial range, you can also specify
    a special entry ``"cname"``, with which parameters can be coupled to each other.
    For example, if two variables both have ``"cname": "SomeName"``, then only one
    of them will be changed, and the other will always have the same value. This can
    come in handy if you'd like to optimize two models, say one for the
    visible matter and one for dark matter, and you'd like both to be centered on
    the exact same position even though this position might need to be determined.
    Perhaps even the orientation could be forced to be the same this way.

    In such a case, it's even possible to perform some operations on the variables:
    if ``"cname": "SomeName ; x*2"`` is specified, then this parameter's value will
    be set to twice the reference value. There should always be one without any extra
    operations, so that a reference value can be determined.

    If a dictionary is used, you can also specify some ``"tag"`` value. This will be
    interpreted as a string, but will mostly be ignored. You can use this to name
    some variables yourself, perhaps to make post-processing lens inversion results
    easier to interpret.
    """
    # TODO: also add { "prior": { "type": "gaussian", "params": [ someMu, someSigma ]}}
    #       to explanation; TODO: is now translated to "neglogprob"

    # TODO:
    #     "neglogprob": [
    #        { "args": [ "var1", "var2" ],
    #          "clcode": "float functionname(float a, float b) { ... ;  return ... }" },
    #        { "args": [ "var1", (123, "var_for_unit") ],
    #          "clcall": "functionname",
    #     ]

    # TODO: "varname"
    # TODO: { "fixed": 123 } for fixed value also possible, allows use of "varname"

    inf = _createTemplateLens(parametricLens, Dd)

    def getNameForOffset(off):
        for x in inf["paramoffsets"]:
            if x["offset"] == off:
                return x["name"]
        return "unknown"

    initMin, initMax = _createInitialMinMaxParameters(inf, defaultFraction, allowSameInitialMinMax)
    hardMin, hardMax = _createHardMinMaxParameters(inf)

    numParams = initMin.shape[0]
    for i in range(numParams):
        if hardMin[i] > initMin[i]:
            if clampToHardLimits:
                initMin[i] = hardMin[i]
            else:
                n = getNameForOffset(i)
                raise ParametricDescriptionException(f"Parameter '{n}' (at offset {i}) has initial value that can be smaller than hard lower limit")
        if hardMax[i] < initMax[i]:
            if clampToHardLimits:
                initMax[i] = hardMax[i]
            else:
                n = getNameForOffset(i)
                raise ParametricDescriptionException(f"Parameter '{n}' (at offset {i}) has initial value that can be larger than hard upper limit")
        if initMax[i] < initMin[i]:
            n = getNameForOffset(i)
            raise ParametricDescriptionException(f"Parameter '{n}' (at offset {i}) has invalid initial range (possibly after clamping): initMin = {initMin[i]}, initMax = {initMax[i]}")

    paramRanges = [ ]
    for offInf in inf["paramoffsets"]:
        offset = offInf["offset"]
        hardMinValue = hardMin[offset]
        hardMaxValue = hardMax[offset]
        initMinValue = initMin[offset]
        initMaxValue = initMax[offset]
        
        assert not (offset in paramRanges), f"Internal error: offset {offset} already set"

        d = { "initialrange": [ initMinValue, initMaxValue ],
                              "hardlimits": [ hardMinValue, hardMaxValue ],
                              "name": offInf["name"],
                              "scalefactor": offInf["scalefactor"],
                              "offset": offset }

        if "cname" in offInf:
            d["cname"] = offInf["cname"]
            d["ccode"] = offInf["ccode"]
            d["cname_orig"] = offInf["cname_orig"]
        if "tag" in offInf:
            d["tag"] = offInf["tag"]
        if "varname" in offInf:
            d["varname"] = offInf["varname"]
        if "prior" in offInf:
            d["prior"] = offInf["prior"]
        paramRanges.append(d)

    paramRanges = sorted(paramRanges, key = lambda x: x["offset"])

    # Do a final consistency check
    for p in paramRanges:
        n = p["name"]
        hardMin, hardMax = p["hardlimits"]
        initMin, initMax = p["initialrange"]
        if initMin == initMax and not allowSameInitialMinMax:
            raise ParametricDescriptionException(f"Consistency error: Parameter '{n}' has emty initial range ({initMin})")

        if not (hardMin <= initMin <= initMax <= hardMax):
            raise ParametricDescriptionException(f"Consistency error: Parameter '{n}' fails range check: hardMin = {hardMin}, initMin = {initMin}, initMax = {initMax}, hardMax = {hardMax}")

    #print("variablefloatparams")
    #pprint.pprint(paramRanges)

    varNameOffsets = { }
    for varNameDict in [ inf["fixedparamoffsets"], paramRanges ]:
        for x in varNameDict:
            if not "varname" in x:
                continue
            varName = x["varname"]
            if varName in varNameOffsets:
                raise ParametricDescriptionException(f"Variable name '{varName}' is used more than once")
            varNameOffsets[varName] = { "offset": x["offset"], "scalefactor": x["scalefactor"] }

    #print("varNameOffsets")
    #pprint.pprint(varNameOffsets)

    clCode = _createOpenClCode(inf["allclcodeinfo"], varNameOffsets, paramRanges)
    
    ret = {
        "templatelens": inf["templatelens"],
        "floatparams": inf["floatparams"],
        "scales": inf["scales"],
        "variablefloatparams": paramRanges,
        "neglogprobcode": clCode
    }
    return ret

def _getFunctionNameAndArguments(code):
    origCode = code
    code = code.replace("\r\n", "\n")
    code = code.replace("\n", " ")
    code = code.strip()
    idx = code.find(" ")
    if idx < 0:
        raise ParametricDescriptionException(f"Can't analyze OpenCL function: {code}")

    retType = code[:idx]
    if retType != "float":
        raise ParametricDescriptionException(f"Expecting return type 'float' in OpenCL code, but is '{retType}' in {code}")

    argStart = code.find("(", idx+1)
    if argStart < 0:
        raise ParametricDescriptionException(f"Can't find start of arguments in OpenCL code {code}")
    funcName = code[idx+1:argStart].strip()
    if not funcName:
        raise ParametricDescriptionException(f"Empty function name in OpenCL code {code}")

    argEnd = code.find(")", argStart+1)
    if argEnd < 0:
        raise ParametricDescriptionException(f"Can't find end of arguments in OpenCL code {code}")
    argString  = code[argStart+1:argEnd].strip()

    if not argString:
        args = []
    else:
        args = argString.split(",")
        args = [ a.strip() for a in args ]
        for a in args:
            if not a:
                raise ParametricDescriptionException(f"Empty argument in OpenCL code {code}")

    return funcName, args

def _createOpenClCode(allClCodeInfo, varNameOffsets, paramRanges):

    functionOrder = [ ]
    uniqueFunctionCode = { }
    negLogProbCalls = [ ]

    # Introduce some defaults
    for knownFunctionsCodes in [
                "float neg_log_gaussian(float x, float mu, float sigma) { float fracdiff = (x-mu)/sigma; return 0.5*fracdiff*fracdiff; }"
            ]:

        funcName, funcArgs = _getFunctionNameAndArguments(knownFunctionsCodes)
        uniqueFunctionCode[funcName] = { "args": funcArgs, "code": knownFunctionsCodes }
        functionOrder.append(funcName)

    # Analyze the explicitly specified OpenCL code
    for inf in allClCodeInfo:

        if "clcode" in inf:
            code = inf["clcode"]
            funcName, funcArgs = _getFunctionNameAndArguments(inf["clcode"])

            if funcName in uniqueFunctionCode:
                # If the same function name is used, and it's actually the same
                # function, don't complain, just ignore. This allows us to introduce
                # a helper function multiple times.
                if uniqueFunctionCode[funcName]["code"] != code:
                    raise ParametricDescriptionException(f"Multiple function definitions for function '{funcName}'")
            else:
                uniqueFunctionCode[funcName] = { "args": funcArgs, "code": code }
                functionOrder.append(funcName)
        else:
            funcName = inf["clcall"]
            if not funcName in uniqueFunctionCode:
                raise ParametricDescriptionException(f"Performing a 'clcall' of function '{funcName}', but this has not been defined")
            funcArgs = uniqueFunctionCode[funcName]["args"]

        if "args" in inf:
            # If this is present, the function should be called with the same
            # number of arguments as 'funcArgs'
            # If not present, it's just a helper function
            argVarNames = inf["args"]
            if len(argVarNames) != len(funcArgs):
                raise ParametricDescriptionException(f"When calling OpenCL function {funcName}, {len(argVarNames)} variables are used, but expecting {len(funcArgs)}")

            offsetsAndValues = []
            for vn in argVarNames:

                if type(vn) == str:
                    if not vn in varNameOffsets:
                        raise ParametricDescriptionException(f"Variable name {vn} in call to OpenCL function {funcName} is not known")

                    offsetsAndValues.append("pFloatParams[{}] /* {} */".format(varNameOffsets[vn]["offset"], vn))
                elif type(vn) == tuple or type(vn) == list:
                    value, scaleVar = vn
                    if not scaleVar in varNameOffsets:
                        raise ParametricDescriptionException(f"Variable name {vn} (for scale factor) in call to OpenCL function {funcName} is not known")

                    value /= varNameOffsets[scaleVar]["scalefactor"]
                    offsetsAndValues.append(json.dumps(value))
                else: # assume it's a float or int value for example
                    value = vn
                    offsetsAndValues.append(json.dumps(value))


            callLine = f"    value += {funcName}(" + ", ".join(offsetsAndValues) + ");";
            negLogProbCalls.append(callLine)

    # Convert the parameter priors into OpenCL calls
    for inf in paramRanges:
        if not "prior" in inf:
            continue
        prior = inf["prior"]
        offset = inf["offset"]
        name = inf["varname"] if "varname" in inf else inf["name"]
        callStr = _toScaledClCall(inf["prior"], f"pFloatParams[{offset}]", inf["scalefactor"])
        callStr += f" // Prior for '{name}', based on: {prior}"

        callLine = f"    value += " + callStr
        negLogProbCalls.append(callLine)

    if not negLogProbCalls:
        return None

    totalCode = ""
    for fn in functionOrder:
        totalCode += uniqueFunctionCode[fn]["code"]
        totalCode += "\n"

    totalCode += """
float calculateTotalNegLogProbContributions(const float *pFloatParams)
{
    float value = 0;
"""
    totalCode += "\n".join(negLogProbCalls)
    totalCode += """
    return value;
}
"""
    return totalCode

def _getUnitlessValue(x):
    return f"{x:.10g}"

def _getUnitValue(x, unitStr):
    unit = eval(unitStr)
    y = x/unit
    v = _getUnitlessValue(y)
    return f"{v}*{unitStr}"

def _convertedValueToString(value, stringConverter):
    if type(value) == list or type(value) == tuple:
        if len(value) == 1:
            return "[" + stringConverter(value[0]) + "]"
        if len(value) == 2:
            return "[" + stringConverter(value[0]) + ", " + _getUnitlessValue(value[1]) + "]"
        if len(value) == 3:
            return "[" + stringConverter(value[0]) + ", " + _getUnitlessValue(value[1]) + ", " + _getUnitlessValue(value[2]) + "]"
        raise ParametricDescriptionException(f"List or tuple should have length 1, 2 or 3, but is {len(value)}")
    
    if type(value) == dict:
        d = "{ "
        for k in value:
            if k == "cname" or k == "tag" or k == "varname":
                cn = str(value[k])
                d += f'"{k}": ' + json.dumps(cn) + ", "
                continue
            if k == "prior":
                d += '"prior": ' + _formatPrior(value[k], stringConverter) + ', '
                continue

            if not k in [ "initmin", "initmax", "hardmin", "hardmax" ]:
                raise ParametricDescriptionException(f'Unexpected key {k}, expecting "initmin", "initmax", "hardmin", "hardmax"')
            d += f'"{k}": ' + stringConverter(value[k]) + ", "
        d += "}"
        return d

    return stringConverter(value)
            
def _analyzePlummerLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):
    params = lens.getLensParameters()
    massStr = _convertedValueToString(convertValueFunction(params["mass"], ["PlummerLens"], "mass", "mass_scaled", params), lambda x: _getUnitValue(x, massUnitString))
    widthStr = _convertedValueToString(convertValueFunction(params["width"], ["PlummerLens"], "width", "width_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "mass": {massStr},',
        f'    "width": {widthStr},',
        '}'
    ]

def _analyzeNSIELens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):
    params = lens.getLensParameters()    
    sigmaStr = _convertedValueToString(convertValueFunction(params["velocityDispersion"], ["NSIELens"], "velocityDispersion", "sigma_scaled", params), _getUnitlessValue)
    ellStr = _convertedValueToString(convertValueFunction(params["ellipticity"], ["NSIELens"], "ellipticity", "ellipticity", params), _getUnitlessValue)
    coreStr = _convertedValueToString(convertValueFunction(params["coreRadius"], ["NSIELens"], "coreRadius", "core_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "velocityDispersion": {sigmaStr},',
        f'    "ellipticity": {ellStr},',
        f'    "coreRadius": {coreStr},',
        '}'
    ]

def _analyzeMultiplePlummerLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):
    paramLines = ['[']
    for i,params in enumerate(lens.getLensParameters()):
        massStr = _convertedValueToString(convertValueFunction(params["mass"], ["MultiplePlummerLens"], f"mass_{i}", f"mass_{i}_scaled", params), lambda x: _getUnitValue(x, massUnitString))
        widthStr = _convertedValueToString(convertValueFunction(params["width"], ["MultiplePlummerLens"], f"width_{i}", f"width_{i}_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
        xStr = _convertedValueToString(convertValueFunction(params["x"], ["MultiplePlummerLens"], f"x_{i}", f"x_{i}_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
        yStr = _convertedValueToString(convertValueFunction(params["y"], ["MultiplePlummerLens"], f"y_{i}", f"y_{i}_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
        subParams = [
            '{',
            f'    "mass": {massStr},',
            f'    "width": {widthStr},',
            f'    "x": {xStr},',
            f'    "y": {yStr},',
            '},' ]
        for p in subParams:
            paramLines.append('    ' + p)

    paramLines.append(']')
    return paramLines

def _analyzeCompositeLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):
    paramLines = ['[']
    for i,params in enumerate(lens.getLensParameters()):
        xStr = _convertedValueToString(convertValueFunction(params["x"], ["CompositeLens"], f"x_{i}", f"x_{i}_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
        yStr = _convertedValueToString(convertValueFunction(params["y"], ["CompositeLens"], f"y_{i}", f"y_{i}_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
        angleStr = _convertedValueToString(convertValueFunction(params["angle"], ["CompositeLens"], f"angle_{i}", f"angle_{i}", params), _getUnitlessValue)
        factorStr = _convertedValueToString(convertValueFunction(params["factor"], ["CompositeLens"], f"factor_{i}", f"factor_{i}", params), _getUnitlessValue)
        subParams = [
            '{',
            f'    "factor": {factorStr},',
            f'    "x": {xStr},',
            f'    "y": {yStr},',
            f'    "angle": {angleStr},' ]

        def cvfWrapper(x, lensName, paramName, uniqueParamName, fullParams):
            return convertValueFunction(x, [ "CompositeLens", f"lens_{i}" ] + lensName, paramName, f"lens_{i}," + uniqueParamName, fullParams)

        subLensLines = createParametricDescription(params["lens"], massUnitString, angularUnitString, False, cvfWrapper, objectStore, objectStoreName)
        subParams.append('    "lens": ' + subLensLines[0])
        for sl in subLensLines[1:]:
            subParams.append('    ' + sl)
        subParams.append('},')

        for p in subParams:
            paramLines.append('    ' + p)

    paramLines.append(']')
    return paramLines

def _analyzeSISLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):
    params = lens.getLensParameters()
    sigmaStr = _convertedValueToString(convertValueFunction(params["velocityDispersion"], ["SISLens"], "velocityDispersion", "sigma_scaled", params), _getUnitlessValue)
    return [
        '{',
        f'    "velocityDispersion": {sigmaStr},',
        '}'
    ]

def _analyzeMassSheetLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):
    params = lens.getLensParameters()
    densStr = _convertedValueToString(convertValueFunction(params["density"], ["MassSheetLens"], "density", "density_scaled", params), _getUnitlessValue)
    return [
        '{',
        f'    "density": {densStr},',
        '}'
    ]

def _analyzeDeflectionGridLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):
    # Nothing can be changed about this lens
    if objectStore is None or objectStoreName is None:
        raise ParametricDescriptionException("To be able to store the deflection angles for a DeflectionGridLens, the objectStore and objectStoreName must be set")

    params = lens.getLensParameters()
    bl = params["bottomleft"]
    tr = params["topright"]
    bl_x = _getUnitValue(bl[0], angularUnitString)
    bl_y = _getUnitValue(bl[1], angularUnitString)
    tr_x = _getUnitValue(tr[0], angularUnitString)
    tr_y = _getUnitValue(tr[1], angularUnitString)

    objectUuid = str(uuid.uuid4())
    objectStore[objectUuid] = params["angles"].copy()

    return [
        '{',
        f'    "bottomleft": [ {bl_x}, {bl_y} ],',
        f'    "topright": [ {tr_x}, {tr_y} ],',
        f'    "angles": {objectStoreName}["{objectUuid}"],',
        '}'
    ]

def _analyzePIMDLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):

    params = lens.getLensParameters()
    sigma0Str = _convertedValueToString(convertValueFunction(params["centraldensity"], ["PIMDLens"], "centraldensity", "centraldensity_scaled", params), _getUnitlessValue)
    coreRadStr = _convertedValueToString(convertValueFunction(params["coreradius"], ["PIMDLens"], "coreradius", "coreradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    scaleRadStr = _convertedValueToString(convertValueFunction(params["scaleradius"], ["PIMDLens"], "scaleradius", "scaleradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "centraldensity": {sigma0Str},',
        f'    "coreradius": {coreRadStr},',
        f'    "scaleradius": {scaleRadStr},',
        '}'
    ]

def _analyzePIEMDLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):

    params = lens.getLensParameters()
    sigma0Str = _convertedValueToString(convertValueFunction(params["centraldensity"], ["PIEMDLens"], "centraldensity", "centraldensity_scaled", params), _getUnitlessValue)
    epsStr = _convertedValueToString(convertValueFunction(params["epsilon"], ["PIEMDLens"], "epsilon", "epsilon", params), _getUnitlessValue)
    coreRadStr = _convertedValueToString(convertValueFunction(params["coreradius"], ["PIEMDLens"], "coreradius", "coreradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    scaleRadStr = _convertedValueToString(convertValueFunction(params["scaleradius"], ["PIEMDLens"], "scaleradius", "scaleradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "centraldensity": {sigma0Str},',
        f'    "epsilon": {epsStr},',
        f'    "coreradius": {coreRadStr},',
        f'    "scaleradius": {scaleRadStr},',
        '}'
    ]

def _analyzeLTPIEMDLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):

    params = lens.getLensParameters()
    velDispStr = _convertedValueToString(convertValueFunction(params["velocitydispersion"], ["LTPIEMDLens"], "velocitydispersion", "velocitydispersion_scaled", params), _getUnitlessValue)
    ellStr = _convertedValueToString(convertValueFunction(params["ellipticity"], ["LTPIEMDLens"], "ellipticity", "ellipticity", params), _getUnitlessValue)
    coreRadStr = _convertedValueToString(convertValueFunction(params["coreradius"], ["LTPIEMDLens"], "coreradius", "coreradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    scaleRadStr = _convertedValueToString(convertValueFunction(params["scaleradius"], ["LTPIEMDLens"], "scaleradius", "scaleradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "velocitydispersion": {velDispStr},',
        f'    "ellipticity": {ellStr},',
        f'    "coreradius": {coreRadStr},',
        f'    "scaleradius": {scaleRadStr},',
        '}'
    ]

def _analyzeLTPIMDLens(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName):

    params = lens.getLensParameters()
    velDispStr = _convertedValueToString(convertValueFunction(params["velocitydispersion"], ["LTPIMDLens"], "velocitydispersion", "velocitydispersion_scaled", params), _getUnitlessValue)
    coreRadStr = _convertedValueToString(convertValueFunction(params["coreradius"], ["LTPIMDLens"], "coreradius", "coreradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    scaleRadStr = _convertedValueToString(convertValueFunction(params["scaleradius"], ["LTPIMDLens"], "scaleradius", "scaleradius_scaled", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "velocitydispersion": {velDispStr},',
        f'    "coreradius": {coreRadStr},',
        f'    "scaleradius": {scaleRadStr},',
        '}'
    ]

# TODO: convert prior parameters to units
def createParametricDescription(lens, massUnitString = "MASS_SUN", angularUnitString = "ANGLE_ARCSEC",
                                asString = True, convertValueFunction = None,
                                objectStore = None, objectStoreName = None):
    """Create a basic representation of a parametric lens model, based on the
    :class:`lens model<grale.lenses.GravitationalLens>` in `lens`. The result is
    a string which represents python code and can be saved to a file for further
    editing. By default this result has no parameters that can change and will need to be
    adjusted, but this behaviour can be changed using the `convertValueFunction`
    parameter (see below).
    
    To make the description more readable, values for masses will be represented
    as a value times the mass unit, which needs to be represented as string in
    `massUnitString`. Similarly, angular values will be represented using 
    `angularUnitString`.

    By default a single large string is returned, if a list of separate lines is
    more covenient, the `asString` parameter can be set to ``False``.

    The `convertValueFunction` is a callback that can be used to convert a fixed
    value to a variable one. The parameters are:
     
        - `value`: the value itself
        - `lensName`: this is a list, of which the last entry specifies the name
          of the current (sub) model. This list could simply be ``[ "NSIELens" ]`` for
          a simple lens, or e.g. ``[ "CompositeLens", "lens_0", "NSIELens" ]`` for a
          more complex model.
        - `paramName`: the name of the parameter, for example ``"mass"``
        - `uniqueParamName`: internal name of the parameter, which is unique. This is
          similar to the parameter names in ``variablefloatparams`` from the
          :func:`analyzeParametricLensDescription` function.
        - `fullParams`: the full lens model parameters of the current (sub) model. This
          could be helpful if you need more information about which submodel is being
          considered.
    
    The return value can be:
     
        - a value: in that case this is a fixed parameter
        - a list ``[ value, fraction ]`` to initialize the parameter to vary between
          bounds specified by the fraction, or just ``[ value ]`` if the default
          fraction is to be used when calling :func:`analyzeParametricLensDescription`.
          If a third fraction is specified, e.g. ``[ value, fraction, hardfraction ]``,
          then that will be used for the hard bounds.
        - a dictionary containing entries for ``"initmin"`` and ``"initmax"``, and
          optionally ``"hardmin"`` and ``"hardmax"``. Entries for ``"cname"`` and
          ``"tag"`` are also allowed - see :func:`analyzeParametricLensDescription`
          for their meaning.

    For some lenses (for now only a :class:`DeflectionGridLens <grale.lenses.DeflectionGridLens>`)
    it may be nessary to store a large amount of data somewhere (e.g. the deflection 
    angles on a grid). In this case you can specify a dictionary for `objectStore` where
    this data will be stored with a unique key. In the final output of this function (which
    is a string), the name `objectStoreName` will be used as the name of this dictionary.
    """

    if convertValueFunction is None:
        convertValueFunction = lambda value, lensName, paramName, uniqueParamName, allParams : value
    
    if not type(lens) in _supportedLensTypesByClass:
        raise ParametricDescriptionException(f"Can't create parametric description for lens type {type(lens)}")
    
    info = _supportedLensTypesByClass[type(lens)]
    name = info["name"]
    if not "analysis" in info:
        raise ParametricDescriptionException(f"No lens analysis function for {name}")
    
    analyzer = info["analysis"]
    paramLines = analyzer(lens, massUnitString, angularUnitString, convertValueFunction, objectStore, objectStoreName)
    descLines = [ 
        '{',
        f'    "type": "{name}",',
        f'    "params": ' + paramLines[0] ]
    for p in paramLines[1:]:
        descLines.append('    ' + p)
    descLines.append('}')
    
    if not asString:
        return descLines
    return "\n".join(descLines)

_supportedLensTypes = {
    "PlummerLens": { "handler": _processPlummerLens, "lens": lenses.PlummerLens, "analysis": _analyzePlummerLens },
    "MultiplePlummerLens": { "handler": _processMultiplePlummerLens, "lens": lenses.MultiplePlummerLens, "analysis": _analyzeMultiplePlummerLens },
    "NSIELens": { "handler": _processNSIELens, "lens": lenses.NSIELens, "analysis": _analyzeNSIELens },
    "CompositeLens": { "handler": _processCompositeLens, "lens": lenses.CompositeLens, "analysis": _analyzeCompositeLens },
    "SISLens": { "handler": _processSISLens, "lens": lenses.SISLens, "analysis": _analyzeSISLens },
    "MassSheetLens": { "handler": _processMassSheetLens, "lens": lenses.MassSheetLens, "analysis": _analyzeMassSheetLens },
    "DeflectionGridLens": { "handler": _processDeflectionGridLens, "lens": lenses.DeflectionGridLens, "analysis": _analyzeDeflectionGridLens },
    "PIMDLens" : { "handler": _processPIMDLens, "lens": lenses.PIMDLens, "analysis": _analyzePIMDLens },
    "PIEMDLens" : { "handler": _processPIEMDLens, "lens": lenses.PIEMDLens, "analysis": _analyzePIEMDLens },
    "LTPIMDLens" : { "handler": _processLTPIMDLens, "lens": lenses.LTPIMDLens, "analysis": _analyzeLTPIMDLens },
    "LTPIEMDLens" : { "handler": _processLTPIEMDLens, "lens": lenses.LTPIEMDLens, "analysis": _analyzeLTPIEMDLens },
}

_supportedLensTypesByClass = { _supportedLensTypes[name]["lens"]: { "name": name, **_supportedLensTypes[name] } for name in _supportedLensTypes }

def getSupportedLensTypes():
    """List which gravitational lens types are recognized in the parametric description."""
    return [ (x, _supportedLensTypes[x]["lens"]) for x in _supportedLensTypes ]

def _refineValue(v1, v2, fraction, lensType, paramName):
    
    v1min, v1max = None, None
    try:
        v1min, v1max = v2*(1-fraction), v2*(1+fraction)
    except:
        # Probably not a number, try a function call
        pass

    if v1min is None or v1max is None:
        v1min, v1max = fraction(v2, lensType, paramName)

    if v1min > v1max:
        v1min, v1max = v1max, v1min

    if type(v1) == list or type(v1) == tuple:
        return { "initmin": v1min, "initmax": v1max }

    if type(v1) != dict: # just a value, keep it
        return v1

    # It's a dict, but still may represent a fixed value
    if "fixed" in v1:
        return v1.copy()

    v1 = v1.copy()
    v1["initmin"] = v1min
    v1["initmax"] = v1max
    return v1

def refineParametricDescription(initialParamDesc, lens, fraction):
    """This is a helper function to adjust an earlier parametric lens description
    based on an estimate of the solution. The envisioned usage is to first do a
    general parametric lens inversion based on `initialParamDesc`, which leads
    to a lens model `lens`. Then you adjust the initial parameter ranges to a narrow
    range around the values from this model, so that a following MCMC exploration can 
    look in the most interesting region of the parameter space.

    Arguments:

     - `initialParamDesc`: the parametric model description (see e.g.
       :func:`analyzeParametricLensDescription`) that is used to obtain an initial
       lens model, which is then passed as the following argument.
     - `lens`: this should be the lens model that results from the parametric optimization
       that was performed with the previous lens description.
       The initial bounds on the parameters that are allowed to change are set to a
       small range around the corresponding value from this model.
     - `fraction`: for the parameters that are not fixed, the initial values will be
       based on this fraction. For example, if this is 0.05, or 5%, the initial value
       range will be set to 95% and 105% of the value of the lens model. Alternatively,
       this can also be a callback function that should return the new initial value
       range. It is called with three arguments, the value itself, the type of the lens
       and the name of the parameter.
    """

    paramDesc = initialParamDesc.copy()
    refineDesc = lens if type(lens) == dict else eval(createParametricDescription(lens)) # TODO: other parameters?
    
    t1, t2 = paramDesc["type"], refineDesc["type"]
    if t1 != t2:
        raise Exception(f'Lens type mismatch: {t1} != {t2}')

    p1, p2 = paramDesc["params"], refineDesc["params"]
    if t1 == "CompositeLens":
        assert type(p1) == list and type(p2) == list, "Parameters of CompositeLens should be lists"
        assert len(p1) == len(p2), "Different number of sublenses in CompositeLens"

        for subP1, subP2 in zip(p1, p2):
            for prop in [ "x", "y", "angle", "factor" ]:
                subP1[prop] = _refineValue(subP1[prop], subP2[prop], fraction, t1, prop)

            subLens1, subLens2 = subP1["lens"], subP2["lens"]
            subLens1 = refineParametricDescription(subLens1, subLens2, fraction)
            subP1["lens"] = subLens1
    else:
        assert type(p1) == dict and type(p2) == dict, "Parameters should be contained in dictionaries"
        for k in p1:
            p1[k] = _refineValue(p1[k], p2[k], fraction, t1, k)

    return paramDesc


def refineParametricDescription_old(varFloatParams, lens, fraction, evaluate = True):
    """This is a helper function to adjust an earlier parametric lens description
    based on an estimate of the solution. The envisioned usage is to first do a
    general parametric lens inversion, which provides `varFloatParams` and leads
    to a lens model `lens`, then adjust the initial parameter ranges to a narrow
    range around the values from this model, so that a following MCMC exploration can 
    look in the most interesting region of the parameter space.
    
    Arguments:
     
     - `varFloatParams`: this should be the description of the variable parameters and
       their bounds from an earlier analysis of a parametric lens description. If the
       first call was made using :func:`analyzeParametricLensDescription`, then the
       returned dictionary contains this input using the ``"variablefloatparams"`` key.
       If the :func:`InversionWorkSpace.invertParametric <grale.inversion.InversionWorkSpace.invertParametric>`
       or :func:`invertParametric <grale.inversion.invertParametric>` calls were used,
       they provide this input as the last entry in the returned tuple.
     - `lens`: this should be the lens model that results from a parametric optimization,
       that was performed with a lens description that led to the `varFloatParams` input.
       The initial bounds on the parameters that are allowed to change are set to a
       small range around the corresponding value from this model.
     - `fraction`: for the parameters that are not fixed, the initial values will be
       based on this fraction. For example, if this is 0.05, or 5%, the initial value
       range will be set to 95% and 105% of the value of the lens model.
     - `evaluate`: by default, the new parametric description will be returned as a Python
       dictionary. If set to `False`, a string will be returned instead.
    """

    paramInfo = { varFloatParams[i]["name"]:varFloatParams[i] for i in range(len(varFloatParams)) }

    def cvf(value, lensName, paramName, uniqueParamName, lensParams):
        if not uniqueParamName in paramInfo:
            return value

        prevParams = paramInfo[uniqueParamName]
        newInitMin = value*(1-fraction)
        newInitMax = value*(1+fraction)
        if newInitMax < newInitMin:
            newInitMin, newInitMax = newInitMax, newInitMin

        scaleFactor = prevParams["scalefactor"]

        newDict = { "initmin": newInitMin, "initmax": newInitMax, 
                 "hardmin": prevParams["hardlimits"][0] * scaleFactor,
                 "hardmax": prevParams["hardlimits"][1] * scaleFactor }

        if "tag" in prevParams:
            newDict["tag"] = prevParams["tag"]
        if "varname" in prevParams:
            newDict["varname"] = prevParams["varname"]
        if "cname_orig" in prevParams:
            newDict["cname"] = prevParams["cname_orig"]
        if "prior" in prevParams:
            newDict["prior"] = prevParams["prior"]

        return newDict

    newDesc = createParametricDescription(lens, convertValueFunction=cvf)
    #print("REFINED DESCRIPTION:")
    #print(newDesc)
    if not evaluate:
        return newDesc
    
    inf = float("inf")
    newDesc = eval(newDesc)
    return newDesc

def _checkGaussianPriorParameters(params):
    if len(params) != 2:
        raise ParametricDescriptionException(f"Gaussian prior should have two parameters, got these: {params}")

    mu, sigma = params
    try:
        # Some operation to check that these are numbers
        mu**2, sigma**2
        return
    except:
        pass
    raise ParametricDescriptionException(f"Either mu ({mu}) or sigma ({sigma}) does not appear to be a number")

def _formatGaussionPriorParameters(params, stringConverter):
    return '[' + stringConverter(params[0]) + ',' + stringConverter(params[1]) + ']'

def _scaleGaussianPriorParameters(params, scaleFactor):
    return [ params[0]/scaleFactor, params[1]/scaleFactor ]

def _scaledClGaussianPriorCode(varString, params):
    return "neg_log_gaussian(" + varString + ", " + json.dumps(params[0]) + ", " + json.dumps(params[1]) + "); "

_supportedPriors = { 
    "gaussian": {
        "check": _checkGaussianPriorParameters,
        "format": _formatGaussionPriorParameters,
        "scale": _scaleGaussianPriorParameters,
        "toscaledclcall": _scaledClGaussianPriorCode
        },
}

def _checkPrior(pr):
    if type(pr) != dict or not 'type' in pr or not 'params' in pr:
        raise ParametricDescriptionException("Prior description needs to be a dictionary with 'type' and 'params' keys")
    pr = copy.copy(pr)

    t = pr["type"]
    del pr["type"]
    params = pr["params"]
    del pr["params"]
    if pr:
        raise ParametricDescriptionException(f"Excess information in prior info: {pr}")

    if not t in _supportedPriors:
        raise ParametricDescriptionException(f"Invalid prior name '{t}'")
    
    checkFunction = _supportedPriors[t]["check"]
    checkFunction(params)

def _formatPrior(pr, stringConverter):
    _checkPrior(pr) # Make sure it's valid before we try to format it
    t = pr["type"]
    formatFunction = _supportedPriors[t]["format"]
    return '{ "type": "' + t + '", "params": ' + formatFunction(pr["params"], stringConverter) + '}'

def _toScaledClCall(prior, varString, scaleFactor):
    _checkPrior(prior)
    t = prior["type"]

    scaleFunction = _supportedPriors[t]["scale"]
    adjParams = scaleFunction(prior["params"], scaleFactor)

    codeFunction = _supportedPriors[t]["toscaledclcall"]
    return codeFunction(varString, adjParams)

def getUnitAdjustedPrior(prior, scaleFactor):
    _checkPrior(prior)
    t = prior["type"]
    scaleFunction = _supportedPriors[t]["scale"]
    return { "type": t, "params": scaleFunction(prior["params"], scaleFactor) }

def main2():
    pprint.pprint(getSupportedLensTypes())

def main():
    Dd = 1000*DIST_MPC
    #l = lenses.PlummerLens(Dd, { "mass": 1e13*MASS_SUN, "width": 1*ANGLE_ARCSEC})
    #l = lenses.NSIELens(Dd, { "velocityDispersion": 100000, "ellipticity": 0.8, "coreRadius": 0.5*ANGLE_ARCSEC})
    # l = lenses.MultiplePlummerLens(Dd, [
    #         { "mass": 0.7e15*MASS_SUN, "width":3*ANGLE_ARCSEC, "x": -2*ANGLE_ARCSEC, "y": 1*ANGLE_ARCSEC },
    #         { "mass": 0.4e15*MASS_SUN, "width":4*ANGLE_ARCSEC, "x": 3*ANGLE_ARCSEC, "y": -1*ANGLE_ARCSEC },
    #     ])
    l = lenses.CompositeLens(Dd, [
        {"lens": lenses.PlummerLens(Dd, { "mass": 0.7e15*MASS_SUN, "width":3*ANGLE_ARCSEC}),
        "x": -2*ANGLE_ARCSEC, "y": 1*ANGLE_ARCSEC, "angle": 10, "factor": 1.1},
        {"lens": lenses.PlummerLens(Dd, { "mass": 0.4e15*MASS_SUN, "width":4*ANGLE_ARCSEC}),
        "x": 3*ANGLE_ARCSEC, "y": -1*ANGLE_ARCSEC, "angle":30, "factor": 0.9}
    ])
    #l = lenses.SISLens(Dd, { "velocityDispersion": 400000 })
    
    #Ds = 1
    #Dds = 0.8
    #l = lenses.MassSheetLens(Dd, { "Ds": Ds, "Dds": Dds })

    # l = lenses.GravitationalLens.load("/home/jori/projects/grale2-git/inversion_examples/example2/inv1.lensdata")

    objStr = createParametricDescription(l)
    print(objStr)
    parametricLens = eval(objStr)
    pprint.pprint(parametricLens)

    inf = analyzeParametricLensDescription(parametricLens, Dd, 0.1)
    pprint.pprint(inf)

    pprint.pprint(inf["templatelens"].getLensParameters())

def main0():

    #parametricLens = {
    #    "type": "PlummerLens",
    #    "params": {
    #        # voor vaste waarde
    #        "mass": [ 1e14*MASS_SUN ],
    #        "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
    #    }
    #}

    # parametricLens = {
    #     "type": "MultiplePlummerLens",
    #     "params": [
    #         { "x": [ 0.1*ANGLE_ARCSEC ], "y": 0, "mass": [1e13*MASS_SUN], "width": 3*ANGLE_ARCSEC },
    #         { "x": 0, "y": 0, "mass": 1e13*MASS_SUN, "width": [ 2*ANGLE_ARCSEC ] }
    #     ]
    # }

    parametricLens = {
       "type": "CompositeLens",
       "params": [
            { "x": 0, "y": [ 0.1*ANGLE_ARCSEC ], "factor": 1, "angle": 0,
              "lens": {
                "type": "PlummerLens",
                "params": {
                    "mass": [ 1e14*MASS_SUN ],
                    "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
                }
              }
            },
            { "x": [0.11*ANGLE_ARCSEC], "y": 0, "factor": [ 1 ], "angle": 0,
              "lens": {
                "type": "PlummerLens",
                "params": {
                    "mass": 1e14*MASS_SUN,
                    "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
                }
              }
            },
            { "x": [0.111*ANGLE_ARCSEC], "y": 0, "factor": 1, 
              "angle": { "initmin": -60, "initmax": 30 },
              "lens": {
                  "type": "NSIELens",
                  "params": {
                      "velocityDispersion": [ 400000 ],
                      "ellipticity": [ 0.8 ],
                      "coreRadius": [ 2*ANGLE_ARCSEC ]
                  }
              }
            },
            { "x": 0, "y": 0, "factor": 1, 
              "angle": 0,
              "lens": {
                  "type": "SISLens",
                  "params": {
                      "velocityDispersion": [ 400000 ],
                  }
              }
            },
            {
              "x": 0, "y": 0, "factor": 1, "angle": 0,
              "lens": {
                  "type": "MassSheetLens",
                  "params": { "density": [ 4.44 ] }
              }
            }
        ]
    }

    inf = analyzeParametricLensDescription(parametricLens, 1000*DIST_MPC, 0.1)
    pprint.pprint(inf)
    

if __name__ == "__main__":
    main()
    
