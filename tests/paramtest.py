from grale.all import *
import sys
import pprint
import copy
import math

def getInitialParameterValue(params, paramKey):
    value = params[paramKey]
    if type(value) == dict:
        return (value["initmin"] + value["initmax"])*0.5, False
    
    if type(value) == list or type(value) == tuple:
        if len(value) == 1 or len(value) == 2:
            return value[0], False
        raise Exception("Too many entries for parameter")
    
    return value, True

def getInitialMinOrMaxParameterValue(params, paramKey, fraction, isMin):
    value = params[paramKey]

    key = "initmin" if isMin else "initmax"
    fracSign = -1 if isMin else 1

    if type(value) == dict:
        return value[key], False
    
    if type(value) == list or type(value) == tuple:
        if len(value) == 1:
            fullFrac = 1 + fracSign*abs(fraction)
        elif len(value) == 2:
            fullFrac = 1 + fracSign*abs(value[1])
        else:
            raise Exception("Too many entries for parameter")
        
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

def processPlummerLens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("mass", "mass_scaled"),
                                                                ("width", "width_scaled")],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for PlummerLens")
    
    return newParams, positionNames

def processCompositeLens(lensParams, Dd, getParamValue):
    compParams = []
    positionNames = []

    for idx, subParams in enumerate(lensParams):
        subParams, newParams = _checkParameterValues(subParams, 
                                                     [ ("x", f"x_{idx}_scaled"),
                                                       ("y", f"y_{idx}_scaled"),
                                                       ("factor", f"factor_{idx}"),
                                                       ("angle", f"angle_{idx}") ],
                                                     getParamValue, positionNames)

        l, posNames = createTemplateLens_helper(subParams["lens"], Dd, getParamValue)
        newParams["lens"] = l
        del subParams["lens"]
        for p in posNames:
            positionNames.append(f"lens_{idx}," + p)
        
        if subParams:
            raise Exception("Excess parameters for CompositeLens")
        
        compParams.append(newParams)
        
    return compParams, positionNames

def processNSIELens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("ellipticity", "ellipticity"),
                                                                ("coreRadius", "core_scaled"),
                                                                ("velocityDispersion", "sigma_scaled") ],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for NSIELens")
    
    return newParams, positionNames

def processSISLens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("velocityDispersion", "sigma_scaled") ],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for SISLens")
    
    return newParams, positionNames

def processMassSheetLens(lensParams, Dd, getParamValue):
    positionNames = []
    if "Ds" in lensParams or "Dds" in lensParams:
        raise Exception("For a mass sheet lens, 'Ds' and 'Dds' cannot be used, use 'density' instead")

    lensParams, newParams = _checkParameterValues(lensParams, [ ("density", "density_scaled") ],
                                                  getParamValue, positionNames)    

    if lensParams:
        raise Exception("Excess parameters for MassSheetLens")
    return newParams, positionNames

def processMultiplePlummerLens(lensParams, Dd, getParamValue):
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

def createTemplateLens_helper(parametricLensDescription, Dd, getParamValue):
    lensType = parametricLensDescription["type"]
    lensParams = parametricLensDescription["params"]
    if not lensType in supportedLensTypes:
        raise Exception("Unknown lens type for parametric inversion")
    
    t = supportedLensTypes[lensType]
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

def createTemplateLens(parametricLensDescription, Dd):

    l, paramNames = createTemplateLens_helper(parametricLensDescription, Dd, getInitialParameterValue)
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

def getShouldBeSameArray(templateLensDescription):
    numParams = templateLensDescription["floatparams"].shape[0]
    shouldBeSame = [ True for _ in range(numParams) ]
    for paramInfo in templateLensDescription["paramoffsets"]:
        offset = paramInfo["offset"]
        assert shouldBeSame[offset], "Internal error: changing same floating point value twice"

        shouldBeSame[offset] = False

    return shouldBeSame

def checkShouldBeSame(templateLensDescription, minParams, maxParams):
    shouldBeSame = getShouldBeSameArray(templateLensDescription)

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

def createInitialMinMaxParameters(templateLensDesciption, defaultFraction):
    
    scales = templateLensDesciption["scales"]
    Dd = templateLensDesciption["templatelens"].getLensDistance()
    parametricLensDescription = templateLensDesciption["description"]
    
    l, paramNames = createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: getInitialMinOrMaxParameterValue(params, key, defaultFraction, True))
    _, initMinParams = l.getCLParameters(**scales)

    l, paramNames = createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: getInitialMinOrMaxParameterValue(params, key, defaultFraction, False))
    _, initMaxParams = l.getCLParameters(**scales)

    assert templateLensDesciption["floatparams"].shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes"
    assert initMaxParams.shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes (2)"

    checkShouldBeSame(templateLensDesciption, initMinParams, initMaxParams)

    return initMinParams, initMaxParams

def getHardMinOrMaxParameterValue(params, paramKey, isMin, knownParamNames, hardInfinite):

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
    
def mergeHardMinOrMaxParameterValue(params, key, knownParamNames, paramOffsets):

    # print("knownParamNames")
    # pprint.pprint(knownParamNames)
    # print("paramOffsets")
    # pprint.pprint(paramOffsets)
    startValue, isFixed = getInitialParameterValue(params, key)
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

def createHardMinMaxParameters(templateLensDesciption):

    # We're going to merge this with hard min/max values
    parametricLensDescription = copy.deepcopy(templateLensDesciption["description"])
    Dd = templateLensDesciption["templatelens"].getLensDistance()
    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    paramOffsets = templateLensDesciption["paramoffsets"]
    paramOffsets = { x["name"]:x for x in paramOffsets }
    scales = templateLensDesciption["scales"]

    l, paramNames = createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda params, key: mergeHardMinOrMaxParameterValue(params, key, knownParamNames, paramOffsets))

    # Now parametricLensDescription is modified, to contain
    # "start", "hardmin" and "hardmax" values for all variable
    # entries

    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    hardInfMinNames = []
    l, _ = createTemplateLens_helper(parametricLensDescription, Dd,
                                     lambda params, key: getHardMinOrMaxParameterValue(params, key, True, knownParamNames, hardInfMinNames))

    _, hardMinParams = l.getCLParameters(**scales) # These should all be finite


    knownParamNames = copy.deepcopy(templateLensDesciption["paramnames"])
    hardInfMaxNames = []
    l, _ = createTemplateLens_helper(parametricLensDescription, Dd,
                                     lambda params, key: getHardMinOrMaxParameterValue(params, key, False, knownParamNames, hardInfMaxNames))

    _, hardMaxParams = l.getCLParameters(**scales) # These should all be finite

    def isNumber(num):
        return not (math.isinf(num) or math.isnan(num))

    for i in range(hardMinParams.shape[0]):
        assert isNumber(hardMinParams[i]), "Internal error: expecting all hard min parameters to start as a real number"

    for i in range(hardMaxParams.shape[0]):
        assert isNumber(hardMaxParams[i]), "Internal error: expecting all hard max parameters to start as a real number"
    
    for paramName in hardInfMinNames:
        offset = paramOffsets[paramName]["offset"]
        hardMinParams[offset] = float("-inf")

    for paramName in hardInfMaxNames:
        offset = paramOffsets[paramName]["offset"]
        hardMaxParams[offset] = float("inf")

    assert templateLensDesciption["floatparams"].shape == hardMinParams.shape, "Internal error: mismatch in floating point parameter shapes (3)"
    assert hardMaxParams.shape == hardMinParams.shape, "Internal error: mismatch in floating point parameter shapes (4)"

    checkShouldBeSame(templateLensDesciption, hardMinParams, hardMaxParams)

    return hardMinParams, hardMaxParams

def analyzeParametricLensDescription(parametricLens, Dd, defaultFraction):
    inf = createTemplateLens(parametricLens, Dd)

    initMin, initMax = createInitialMinMaxParameters(inf, defaultFraction)
    hardMin, hardMax = createHardMinMaxParameters(inf)

    numParams = initMin.shape[0]
    for i in range(numParams):
        if hardMin[i] > initMin[i]:
            raise Exception(f"Parameter at offset {i} has initial value that can be smaller than hard lower limit")
        if hardMax[i] < initMax[i]:
            raise Exception(f"Parameter at offset {i} has initial value that can be larger than hard upper limit")

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

def _getUnitValue(x, unitStr):
    unit = eval(unitStr)
    y = x/unit
    return f"{y:.10g}*{unitStr}"

def analyzePlummerLens(lens, massUnitString, angularUnitString):
    params = lens.getLensParameters()
    massStr = _getUnitValue(params["mass"], massUnitString)
    widthStr = _getUnitValue(params["width"], angularUnitString)
    return [
        '{',
        f'    "mass": {massStr},',
        f'    "width": {widthStr}',
        '}'
    ]

def createParametricDescription(lens, massUnitString = "MASS_SUN", angularUnitString = "ANGLE_ARCSEC"):
    
    lookup = { supportedLensTypes[name]["lens"]: { "name": name, **supportedLensTypes[name] } for name in supportedLensTypes }
    if not type(lens) in lookup:
        raise Exception(f"Can't create parametric description for lens type {type(lens)}")
    
    info = lookup[type(lens)]
    if not "analysis" in info:
        raise Exception(f"No lens analysis function for {info['name']}")
    
    analyzer = info["analysis"]
    return analyzer(lens, massUnitString, angularUnitString)

supportedLensTypes = {
    "PlummerLens": { "handler": processPlummerLens, "lens": lenses.PlummerLens, "analysis": analyzePlummerLens },
    "MultiplePlummerLens": { "handler": processMultiplePlummerLens, "lens": lenses.MultiplePlummerLens },
    "NSIELens": { "handler": processNSIELens, "lens": lenses.NSIELens },
    "CompositeLens": { "handler": processCompositeLens, "lens": lenses.CompositeLens },
    "SISLens": { "handler": processSISLens, "lens": lenses.SISLens },
    "MassSheetLens": { "handler": processMassSheetLens, "lens": lenses.MassSheetLens },
}

def main():
    l = lenses.PlummerLens(1000*DIST_MPC, { "mass": 1e13*MASS_SUN, "width": 1*ANGLE_ARCSEC})
    lines = createParametricDescription(l)
    print("\n".join(lines))

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
    