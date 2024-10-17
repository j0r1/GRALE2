"""This module contains tools for parametric inversion: something to analyze
a description of a lens model that can be optimized parametrically, and
a routine to start such a description based on an existing lens model."""

# TODO: separate class for Exception

from .constants import *
from . import lenses
import pprint
import copy
import math

def _getInitialParameterValue(params, paramKey):
    value = params[paramKey]
    if type(value) == dict:
        return (value["initmin"] + value["initmax"])*0.5, False
    
    if type(value) == list or type(value) == tuple:
        if len(value) == 1 or len(value) == 2:
            return value[0], False
        raise Exception("Too many entries for parameter")
    
    return value, True

def _getInitialMinOrMaxParameterValue(params, paramKey, fraction, isMin):
    value = params[paramKey]

    key = "initmin" if isMin else "initmax"
    fracSign = -1 if isMin else 1

    if type(value) == dict:
        return value[key], False
    
    if type(value) == list or type(value) == tuple:

        if len(value) < 1 or len(value) > 2:
            raise Exception("Incorrect number of entries for parameter")

        if value[0] < 0:
            fracSign = -fracSign

        if len(value) == 1:
            fullFrac = 1 + fracSign*abs(fraction)
        else: # 2 params
            fullFrac = 1 + fracSign*abs(value[1])
        
        return value[0]*fullFrac, False
        
    return value, True

def _checkParameterValues(lensParams, paramMapping, getParamValue, positionNames):
    remainingLensParams = lensParams.copy() # shallow copy is enough
    newParams = { }
    
    for x,y in paramMapping:
        
        newParams[x], isFixed = getParamValue(lensParams, x)
        del remainingLensParams[x]
        if not isFixed:
            positionNames.append(y)

    return remainingLensParams, newParams

def _processPlummerLens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("mass", "mass_scaled"),
                                                                ("width", "width_scaled")],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for PlummerLens")
    
    return newParams, positionNames

def _processCompositeLens(lensParams, Dd, getParamValue):
    compParams = []
    positionNames = []

    for idx, subParams in enumerate(lensParams):
        subParams, newParams = _checkParameterValues(subParams, 
                                                     [ ("x", f"x_{idx}_scaled"),
                                                       ("y", f"y_{idx}_scaled"),
                                                       ("factor", f"factor_{idx}"),
                                                       ("angle", f"angle_{idx}") ],
                                                     getParamValue, positionNames)

        l, posNames = _createTemplateLens_helper(subParams["lens"], Dd, getParamValue)
        newParams["lens"] = l
        del subParams["lens"]
        for p in posNames:
            positionNames.append(f"lens_{idx}," + p)
        
        if subParams:
            raise Exception("Excess parameters for CompositeLens")
        
        compParams.append(newParams)
        
    return compParams, positionNames

def _processNSIELens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("ellipticity", "ellipticity"),
                                                                ("coreRadius", "core_scaled"),
                                                                ("velocityDispersion", "sigma_scaled") ],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for NSIELens")
    
    return newParams, positionNames

def _processSISLens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("velocityDispersion", "sigma_scaled") ],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for SISLens")
    
    return newParams, positionNames

def _processMassSheetLens(lensParams, Dd, getParamValue):
    positionNames = []
    if "Ds" in lensParams or "Dds" in lensParams:
        raise Exception("For a mass sheet lens, 'Ds' and 'Dds' cannot be used, use 'density' instead")

    lensParams, newParams = _checkParameterValues(lensParams, [ ("density", "density_scaled") ],
                                                  getParamValue, positionNames)    

    if lensParams:
        raise Exception("Excess parameters for MassSheetLens")
    return newParams, positionNames

def _processMultiplePlummerLens(lensParams, Dd, getParamValue):
    positionNames = []
    newMultiParams = []
    for i,subParams in enumerate(lensParams):
        subParams, newParams = _checkParameterValues(subParams, [ ("x", f"x_{i}_scaled"),
                                                                  ("y", f"y_{i}_scaled"),
                                                                  ("mass", f"mass_{i}_scaled"),
                                                                  ("width", f"width_{i}_scaled"),
                                                                  ],
                                                     getParamValue, positionNames)
        if subParams:
            raise Exception("Excess parameters for MultiplePlummerLens")
        
        newMultiParams.append(newParams)

    return newMultiParams, positionNames

def _createTemplateLens_helper(parametricLensDescription, Dd, getParamValue):
    lensType = parametricLensDescription["type"]
    lensParams = parametricLensDescription["params"]
    if not lensType in _supportedLensTypes:
        raise Exception("Unknown lens type for parametric inversion")
    
    t = _supportedLensTypes[lensType]
    handler, lensClass = t["handler"], t["lens"]

    newParams, positionNames = handler(lensParams, Dd, getParamValue) # Need Dd as a parameter because of possible recursion
    lens = lensClass(Dd, newParams)

    return lens, positionNames        

def _createParamOffsetInfo(l, paramNames):
    adjustableParams = l.getCLAdjustableFloatingPointParameterInfo()
    adjustableParamsDict = { }
    for p in adjustableParams:
        name = p["name"]
        if name in adjustableParamsDict:
            raise Exception(f"Internal error: name {name} already exists in adjustable params dict")
        adjustableParamsDict[name] = p

    paramOffsetInfo = []    
    for n in paramNames:
        if not n in adjustableParamsDict:
            raise Exception(f"Internal error: specified adjustable param {n} does not seem to be valid")
        paramOffsetInfo.append(adjustableParamsDict[n])

    return sorted(paramOffsetInfo, key=lambda x: x["offset"])

def _createTemplateLens(parametricLensDescription, Dd):

    l, paramNames = _createTemplateLens_helper(parametricLensDescription, Dd, _getInitialParameterValue)
    scales = l.getSuggestedScales()
    intParam, floatParams = l.getCLParameters(**scales)
    ret = { "templatelens": l,
            "paramoffsets": _createParamOffsetInfo(l, paramNames),
            "paramnames": paramNames, # in original order
            "scales": scales,
            "floatparams": floatParams,
            "description": copy.deepcopy(parametricLensDescription)
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

def _checkShouldBeSame(templateLensDescription, minParams, maxParams):
    shouldBeSame = _getShouldBeSameArray(templateLensDescription)

    for idx, value in enumerate(shouldBeSame):
        # print(f"Comparing index {idx}: {initMinParams[idx]} vs {initMaxParams[idx]}")
        if value: # everything should be same
            assert minParams[idx] == maxParams[idx], "Internal error: min/max should be same at this position"
            assert minParams[idx] == templateLensDescription["floatparams"][idx], "Internal error: min and template values should match at this position"

        else:
            if minParams[idx] == maxParams[idx]:
                raise Exception(f"Values at floating point offset {idx} should be allowed to change, but initial min/max values are the same (so no variation will be introduced)")

            if minParams[idx] > maxParams[idx]:
                raise Exception(f"Min/max value at floating point offset {idx} should be other way around")

def _createInitialMinMaxParameters(templateLensDesciption, defaultFraction):
    
    scales = templateLensDesciption["scales"]
    Dd = templateLensDesciption["templatelens"].getLensDistance()
    parametricLensDescription = templateLensDesciption["description"]
    
    l, paramNames = _createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: _getInitialMinOrMaxParameterValue(params, key, defaultFraction, True))
    _, initMinParams = l.getCLParameters(**scales)

    l, paramNames = _createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: _getInitialMinOrMaxParameterValue(params, key, defaultFraction, False))
    _, initMaxParams = l.getCLParameters(**scales)

    assert templateLensDesciption["floatparams"].shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes"
    assert initMaxParams.shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes (2)"

    _checkShouldBeSame(templateLensDesciption, initMinParams, initMaxParams)

    return initMinParams, initMaxParams

def _getHardMinOrMaxParameterValue(params, paramKey, isMin, knownParamNames, hardInfinite):

    # We're expecting either a fixed value, or a dict { "start", "hardmin", "hardmax"}
    value = params[paramKey]
    if type(value) != dict:
        return value, True

    paramName = knownParamNames.pop(0)

    key, infValue = ("hardmin", float("-inf")) if isMin else ("hardmax", float("inf"))
    extValue = value[key]
    if math.isinf(extValue):
        paramValue = value["start"] # Don't use inf as a lens parameter
        if extValue != infValue:
            raise Exception(f"Got {extValue} but expecting {infValue} for {paramName}")
        hardInfinite.append(paramName) # Remember this parameter name
    else:
        paramValue = extValue

    return paramValue, False
    
def _mergeHardMinOrMaxParameterValue(params, key, knownParamNames, paramOffsets):

    # print("knownParamNames")
    # pprint.pprint(knownParamNames)
    # print("paramOffsets")
    # pprint.pprint(paramOffsets)
    startValue, isFixed = _getInitialParameterValue(params, key)
    if isFixed:
        return startValue, True

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

    l, paramNames = _createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: _mergeHardMinOrMaxParameterValue(params, key, knownParamNames, paramOffsets))

    # Now parametricLensDescription is modified, to contain
    # "start", "hardmin" and "hardmax" values for all variable
    # entries

    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    hardInfMinNames = []
    l, _ = _createTemplateLens_helper(parametricLensDescription, Dd,
                                     lambda params, key: _getHardMinOrMaxParameterValue(params, key, True, knownParamNames, hardInfMinNames))

    _, hardMinParams = l.getCLParameters(**scales) # These should all be finite


    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    hardInfMaxNames = []
    l, _ = _createTemplateLens_helper(parametricLensDescription, Dd,
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

def analyzeParametricLensDescription(parametricLens, Dd, defaultFraction, clampToHardLimits = False):
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
               # Again a fixed value
               "width": 3*ANGLE_ARCSEC
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
    """
    inf = _createTemplateLens(parametricLens, Dd)

    def getNameForOffset(off):
        for x in inf["paramoffsets"]:
            if x["offset"] == off:
                return x["name"]
        return "unknown"

    initMin, initMax = _createInitialMinMaxParameters(inf, defaultFraction)
    hardMin, hardMax = _createHardMinMaxParameters(inf)

    numParams = initMin.shape[0]
    for i in range(numParams):
        if hardMin[i] > initMin[i]:
            if clampToHardLimits:
                initMin[i] = hardMin[i]
            else:
                n = getNameForOffset(i)
                raise Exception(f"Parameter '{n}' (at offset {i}) has initial value that can be smaller than hard lower limit")
        if hardMax[i] < initMax[i]:
            if clampToHardLimits:
                initMax[i] = hardMax[i]
            else:
                n = getNameForOffset(i)
                raise Exception(f"Parameter '{n}' (at offset {i}) has initial value that can be larger than hard upper limit")

    paramRanges = [ ]
    for offInf in inf["paramoffsets"]:
        offset = offInf["offset"]
        hardMinValue = hardMin[offset]
        hardMaxValue = hardMax[offset]
        initMinValue = initMin[offset]
        initMaxValue = initMax[offset]
        
        assert not (offset in paramRanges), f"Internal error: offset {offset} already set"
        paramRanges.append({ "initialrange": [ initMinValue, initMaxValue ],
                              "hardlimits": [ hardMinValue, hardMaxValue ],
                              "name": offInf["name"],
                              "offset": offset })
    paramRanges = sorted(paramRanges, key = lambda x: x["offset"])
    
    ret = {
        "templatelens": inf["templatelens"],
        "floatparams": inf["floatparams"],
        "scales": inf["scales"],
        "variablefloatparams": paramRanges
    }
    return ret

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
        raise Exception(f"List or tuple should have length 1 or 2, but is {len(value)}")
    
    if type(value) == dict:
        d = "{ "
        for k in value:
            if not k in [ "initmin", "initmax", "hardmin", "hardmax" ]:
                raise Exception(f'Unexpected key {k}, expecting "initmin", "initmax", "hardmin", "hardmax"')
            d += f'"{k}": ' + stringConverter(value[k]) + ", "
        d += "}"
        return d

    return stringConverter(value)
            
def _analyzePlummerLens(lens, massUnitString, angularUnitString, convertValueFunction):
    params = lens.getLensParameters()
    massStr = _convertedValueToString(convertValueFunction(params["mass"], ["PlummerLens"], "mass", params), lambda x: _getUnitValue(x, massUnitString))
    widthStr = _convertedValueToString(convertValueFunction(params["width"], ["PlummerLens"], "width", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "mass": {massStr},',
        f'    "width": {widthStr},',
        '}'
    ]

def _analyzeNSIELens(lens, massUnitString, angularUnitString, convertValueFunction):
    params = lens.getLensParameters()    
    sigmaStr = _convertedValueToString(convertValueFunction(params["velocityDispersion"], ["NSIELens"], "velocityDispersion", params), _getUnitlessValue)
    ellStr = _convertedValueToString(convertValueFunction(params["ellipticity"], ["NSIELens"], "ellipticity", params), _getUnitlessValue)
    coreStr = _convertedValueToString(convertValueFunction(params["coreRadius"], ["NSIELens"], "coreRadius", params), lambda x: _getUnitValue(x, angularUnitString))
    return [
        '{',
        f'    "velocityDispersion": {sigmaStr},',
        f'    "ellipticity": {ellStr},',
        f'    "coreRadius": {coreStr},',
        '}'
    ]

def _analyzeMultiplePlummerLens(lens, massUnitString, angularUnitString, convertValueFunction):
    paramLines = ['[']
    for i,params in enumerate(lens.getLensParameters()):
        massStr = _convertedValueToString(convertValueFunction(params["mass"], ["MultiplePlummerLens"], f"mass_{i}", params), lambda x: _getUnitValue(x, massUnitString))
        widthStr = _convertedValueToString(convertValueFunction(params["width"], ["MultiplePlummerLens"], f"width_{i}", params), lambda x: _getUnitValue(x, angularUnitString))
        xStr = _convertedValueToString(convertValueFunction(params["x"], ["MultiplePlummerLens"], f"x_{i}", params), lambda x: _getUnitValue(x, angularUnitString))
        yStr = _convertedValueToString(convertValueFunction(params["y"], ["MultiplePlummerLens"], f"y_{i}", params), lambda x: _getUnitValue(x, angularUnitString))
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

def _analyzeCompositeLens(lens, massUnitString, angularUnitString, convertValueFunction):
    paramLines = ['[']
    for i,params in enumerate(lens.getLensParameters()):
        xStr = _convertedValueToString(convertValueFunction(params["x"], ["CompositeLens"], f"x_{i}", params), lambda x: _getUnitValue(x, angularUnitString))
        yStr = _convertedValueToString(convertValueFunction(params["y"], ["CompositeLens"], f"y_{i}", params), lambda x: _getUnitValue(x, angularUnitString))
        angleStr = _convertedValueToString(convertValueFunction(params["angle"], ["CompositeLens"], f"angle_{i}", params), _getUnitlessValue)
        factorStr = _convertedValueToString(convertValueFunction(params["factor"], ["CompositeLens"], f"factor_{i}", params), _getUnitlessValue)
        subParams = [
            '{',
            f'    "factor": {factorStr},',
            f'    "x": {xStr},',
            f'    "y": {yStr},',
            f'    "angle": {angleStr},' ]

        def cvfWrapper(x, lensName, paramName, fullParams):
            return convertValueFunction(x, [ "CompositeLens", f"lens_{i}" ] + lensName, paramName, fullParams)

        subLensLines = createParametricDescription(params["lens"], massUnitString, angularUnitString, False, cvfWrapper)
        subParams.append('    "lens": ' + subLensLines[0])
        for sl in subLensLines[1:]:
            subParams.append('    ' + sl)
        subParams.append('},')

        for p in subParams:
            paramLines.append('    ' + p)

    paramLines.append(']')
    return paramLines

def _analyzeSISLens(lens, massUnitString, angularUnitString, convertValueFunction):
    params = lens.getLensParameters()
    sigmaStr = _convertedValueToString(convertValueFunction(params["velocityDispersion"], ["SISLens"], "velocityDispersion", params), _getUnitlessValue)
    return [
        '{',
        f'    "velocityDispersion": {sigmaStr},',
        '}'
    ]

def _analyzeMassSheetLens(lens, massUnitString, angularUnitString, convertValueFunction):
    params = lens.getLensParameters()
    densStr = _convertedValueToString(convertValueFunction(params["density"], ["MassSheetLens"], "density", params), _getUnitlessValue)
    return [
        '{',
        f'    "density": {densStr},',
        '}'
    ]

def createParametricDescription(lens, massUnitString = "MASS_SUN", angularUnitString = "ANGLE_ARCSEC",
                                asString = True, convertValueFunction = None):
    """Create a basic representation of a parametric lens model, based on the
    :class:`lens model<grale.lenses.GravitationalLens>` in `lens`. The result is
    a string which represents python code and can be saved to a file for further
    editing. This result has no parameters that can change, so it will need to be
    adjusted.
    
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
        - `fullParams`: the full lens model parameters of the current (sub) model. This
          could be helpful if you need more information about which submodel is being
          considered.
    
    The return value can be:
     
        - a value: in that case this is a fixed parameter
        - a list ``[ value, fraction ]`` to initialize the parameter to vary between
          bounds specified by the fraction, or just ``[ value ]`` if the default
          fraction is to be used when calling :func:`analyzeParametricLensDescription`
        - a dictionary containing entries for ``"initmin"`` and ``"initmax"``, and
          optionally ``"hardmin"`` and ``"hardmax"``.
    """

    if convertValueFunction is None:
        convertValueFunction = lambda value, lensName, paramName, allParams : value
    
    if not type(lens) in _supportedLensTypesByClass:
        raise Exception(f"Can't create parametric description for lens type {type(lens)}")
    
    info = _supportedLensTypesByClass[type(lens)]
    name = info["name"]
    if not "analysis" in info:
        raise Exception(f"No lens analysis function for {name}")
    
    analyzer = info["analysis"]
    paramLines = analyzer(lens, massUnitString, angularUnitString, convertValueFunction)
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
}

_supportedLensTypesByClass = { _supportedLensTypes[name]["lens"]: { "name": name, **_supportedLensTypes[name] } for name in _supportedLensTypes }

def getSupportedLensTypes():
    """List which gravitational lens types are recognized in the parametric description."""
    return [ (x, _supportedLensTypes[x]["lens"]) for x in _supportedLensTypes ]

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
    
