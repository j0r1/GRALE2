from grale.all import *
import sys
import pprint
import copy

def getInitialParameterValue(value):
    if type(value) == dict:
        return (value["initmin"] + value["initmax"])*0.5, False
    
    if type(value) == list or type(value) == tuple:
        if len(value) == 1 or len(value) == 2:
            return value[0], False
        raise Exception("Too many entries for parameter")
    
    return value, True

def getMinOrMaxParameterValue(value, fraction, keyPrefix, isMin, allowMissingKey):
    key = keyPrefix + "min" if isMin else keyPrefix + "max"
    fracSign = -1 if isMin else 1

    if type(value) == dict:
        if not key in value and allowMissingKey: # for the hard limits
            return None, False
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
    lensParams = lensParams.copy() # shallow copy is enough
    newParams = { }
    
    for x,y in paramMapping:
        
        newParams[x], isFixed = getParamValue(lensParams[x])
        del lensParams[x]
        if not isFixed:
            positionNames.append(y)

    return lensParams, newParams

def processPlummerLens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("mass", "mass_scaled"),
                                                                ("width", "width_scaled")],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for PlummerLens")
    
    return lenses.PlummerLens(Dd, newParams), positionNames

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
        
    return lenses.CompositeLens(Dd, compParams), positionNames

def processNSIELens(lensParams, Dd, getParamValue):
    positionNames = [ ]
    lensParams, newParams = _checkParameterValues(lensParams, [ ("ellipticity", "ellipticity"),
                                                                ("coreRadius", "core_scaled"),
                                                                ("velocityDispersion", "sigma_scaled") ],
                                                  getParamValue, positionNames)    
    if lensParams:
        raise Exception("Excess parameters for NSIELens")
    
    return lenses.NSIELens(Dd, newParams), positionNames

def createTemplateLens_helper(parametricLensDescription, Dd, getParamValue):
    lensType = parametricLensDescription["type"]
    lensParams = parametricLensDescription["params"]
    if lensType == "PlummerLens":
        return processPlummerLens(lensParams, Dd, getParamValue)
    elif lensType == "MultiplePlummerLens":
        pass
    elif lensType == "NSIELens":
        return processNSIELens(lensParams, Dd, getParamValue)
    elif lensType == "CompositeLens":
        return processCompositeLens(lensParams, Dd, getParamValue)
    elif lensType == "SISLens":
        pass
    elif lensType == "MassSheetLens":
        pass
    else:
        raise Exception("Unknown lens type for parametric inversion")

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
            "scales": scales,
            "floatparams": floatParams,
            "description": copy.deepcopy(parametricLensDescription)
          }
    return ret

def createInitialMinMaxParameters(templateLensDesciption, defaultFraction):
    
    scales = templateLensDesciption["scales"]
    Dd = templateLensDesciption["templatelens"].getLensDistance()
    parametricLensDescription = templateLensDesciption["description"]
    
    l, paramNames = createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda value: getMinOrMaxParameterValue(value, defaultFraction, "init", True, False))
    _, initMinParams = l.getCLParameters(**scales)

    l, paramNames = createTemplateLens_helper(parametricLensDescription, Dd,
                                              lambda value: getMinOrMaxParameterValue(value, defaultFraction, "init", False, False))
    _, initMaxParams = l.getCLParameters(**scales)

    assert templateLensDesciption["floatparams"].shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes"
    assert initMaxParams.shape == initMinParams.shape, "Internal error: mismatch in floating point parameter shapes (2)"

    numParams = initMaxParams.shape[0]
    shouldBeSame = [ True for _ in range(numParams) ]
    for paramInfo in templateLensDesciption["paramoffsets"]:
        offset = paramInfo["offset"]
        assert shouldBeSame[offset], "Internal error: changing same floating point value twice"

        shouldBeSame[offset] = False

    # print("paramInfo:")
    # print(paramInfo)
    # print("Min:")
    # print(initMinParams)
    # print("Max:")
    # print(initMaxParams)
    for idx, value in enumerate(shouldBeSame):
        # print(f"Comparing index {idx}: {initMinParams[idx]} vs {initMaxParams[idx]}")
        if value: # everything should be same
            assert initMinParams[idx] == initMaxParams[idx], "Internal error: min/max should be same at this position"
            assert initMinParams[idx] == templateLensDesciption["floatparams"][idx], "Internal error: min and template values should match at this position"

        else:
            if initMinParams[idx] == initMaxParams[idx]:
                raise Exception(f"Values at floating point offset {idx} should be allowed to change, but initial min/max values are the same (so no variation will be introduced)")

    return initMinParams, initMaxParams

def main():

    #parametricLens = {
    #    "type": "PlummerLens",
    #    "params": {
    #        # voor vaste waarde
    #        "mass": [ 1e14*MASS_SUN ],
    #        "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
    #    }
    #}

    parametricLens = {
       "type": "CompositeLens",
       "params": [
            { "x": 0, "y": [ 0.1 ], "factor": 1, "angle": 0,
              "lens": {
                "type": "PlummerLens",
                "params": {
                    "mass": [ 1e14*MASS_SUN ],
                    "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
                }
              }
            },
            { "x": [0.11], "y": 0, "factor": [ 1 ], "angle": 0,
              "lens": {
                "type": "PlummerLens",
                "params": {
                    "mass": 1e14*MASS_SUN,
                    "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
                }
              }
            },
            { "x": [0.111], "y": 0, "factor": 1, 
              "angle": { "initmin": 0, "initmax": 3.14 },
              "lens": {
                  "type": "NSIELens",
                  "params": {
                      "velocityDispersion": [ 400000 ],
                      "ellipticity": [ 0.8 ],
                      "coreRadius": [ 2*ANGLE_ARCSEC ]
                  }
              }
            }
        ]
    }
    
    inf = createTemplateLens(parametricLens, 1000*DIST_MPC)
    print("Template lens info")
    pprint.pprint(inf)
    minMax = createInitialMinMaxParameters(inf, 0.1)
    print("initial min/max info")
    pprint.pprint(minMax)
    

if __name__ == "__main__":
    main()
    





