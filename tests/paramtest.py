from grale.all import *
import sys
import pprint

def getInitialParameterValue(value):
    if type(value) == dict:
        return (value["initmin"] + value["initmax"])*0.5, False
    
    if type(value) == list or type(value) == tuple:
        if len(value) == 1 or len(value) == 2:
            return value[0], False
        raise Exception("Too many entries for parameter")
    
    return value, True

def getMinOrMaxParameterValue(value, fraction, keyPrefix, isMin):
    key = keyPrefix + "min" if isMin else keyPrefix + "max"
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

def createTemplateLens(parametricLensDescription, Dd):

    return createTemplateLens_helper(parametricLensDescription, Dd, getInitialParameterValue)
    
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
            { "x": 0, "y": [ 0 ], "factor": 1, "angle": 0,
              "lens": {
                "type": "PlummerLens",
                "params": {
                    "mass": [ 1e14*MASS_SUN ],
                    "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
                }
              }
            },
            { "x": [0], "y": 0, "factor": [ 1 ], "angle": 0,
              "lens": {
                "type": "PlummerLens",
                "params": {
                    "mass": 1e14*MASS_SUN,
                    "width": { "initmin": 2*ANGLE_ARCSEC, "initmax": 5*ANGLE_ARCSEC }
                }
              }
            },
            { "x": [0], "y": 0, "factor": 1, "angle": [ 0 ],
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
    l, p = createTemplateLens(parametricLens, 1000*DIST_MPC)

    print("Lens")
    print(l)
    print("OpenCL variable param names:")
    pprint.pprint(p)
    print("Adjustable parameters:")
    pprint.pprint(l.getCLAdjustableFloatingPointParameterInfo())

if __name__ == "__main__":
    main()
    





