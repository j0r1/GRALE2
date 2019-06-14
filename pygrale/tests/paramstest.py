from grale.lenses import *
from grale.constants import *
import copy
import pprint
import numpy as np

Dd = 1000*DIST_MPC
AS = ANGLE_ARCSEC

testLenses = [
        {
            "name": "GaussLens",
            "type": GaussLens,
            "params": { "mass": 1e14*MASS_SUN,
                        "width": 10*AS },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "MultiplePlummerLens",
            "type": MultiplePlummerLens,
            "params": [ { "x": -2*AS, "y": 0, "mass": 1e14*MASS_SUN, "width": 10*AS },
                        { "x": 0, "y": 2*AS, "mass": 0.5e14*MASS_SUN, "width": 7*AS } ],
            "direct": True,
            "subpolynomial": False
        },
        { 
            "name": "PlummerLens",
            "type": PlummerLens,
            "params": { "mass": 1e14*MASS_SUN,
                        "width": 10*AS },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "PointmassLens",
            "type": PointmassLens,
            "params": { "mass": 1e14*MASS_SUN },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "SISLens",
            "type": SISLens,
            "params": { "velocityDispersion": 300000 },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "NSIELens",
            "type": NSIELens,
            "params": { "velocityDispersion": 300000, "ellipticity": 0.8, "coreRadius": 1.0*AS },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "NSISLens",
            "type": NSISLens,
            "params": { "velocityDispersion": 300000, "coreRadius": 1.0*AS },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "SIELens",
            "type": SIELens,
            "params": { "velocityDispersion": 300000, "ellipticity": 0.8 },
            "direct": True,
            "subpolynomial": False
        },
        { 
            "name": "SquareLens",
            "type": SquareLens,
            "params": { "mass": 1e14*MASS_SUN,
                        "width": 10*AS },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "MultipleSquareLens",
            "type": MultipleSquareLens,
            "params": [ { "x": -2*AS, "y": 0, "mass": 1e14*MASS_SUN, "width": 10*AS },
                        { "x": 0, "y": 2*AS, "mass": 0.5e14*MASS_SUN, "width": 7*AS } ],
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "MultipleGaussLens",
            "type": MultipleGaussLens,
            "params": [ { "x": -2*AS, "y": 0, "mass": 1e14*MASS_SUN, "width": 10*AS },
                        { "x": 0, "y": 2*AS, "mass": 0.5e14*MASS_SUN, "width": 7*AS } ],
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "MassSheetLens",
            "type": MassSheetLens,
            "params": { "density": 1.1 },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "MassSheetLens",
            "type": MassSheetLens,
            "params": { "Dd": Dd, "Ds": 1, "Dds": 0.8 },
            "direct": False,
            "subpolynomial": False
        },
        {
            "name": "MassDiskLens",
            "type": MassDiskLens,
            "params": { "density": 1.1, "radius": 10*AS },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "MassDiskLens",
            "type": MassDiskLens,
            "params": { "Dd": Dd, "Ds": 1, "Dds": 0.8, "radius": 10*AS },
            "direct": False,
            "subpolynomial": False
        },
        {
            "name": "PolynomialMassProfileLens",
            "type": PolynomialMassProfileLens,
            "params": [
                { "xoffset": 0, "yoffset": 1, "xscale": 5, "yscale": 2, "xend": 3, "coeffs": [ 2, 1, 0.1] },
                { "xoffset": 3, "yoffset": 2, "xscale": 3, "yscale": 4, "xend": 10, "coeffs": [ 1, 2] },
            ],
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "TimeDelayAdjustLens",
            "type": TimeDelayAdjustLens,
            "params": { "z": 1.0, "theta1": 1*AS, "theta2": 5*AS, "dt": 60*60*24 },
            "direct": False,
            "subpolynomial": True
        },
        {
            "name": "ZeroMassLens",
            "type": ZeroMassLens,
            "params": { "density": 1.2, "radius": 3*AS, "zeropoint": 1*AS },
            "direct": False,
            "subpolynomial": True
        },
        {
            "name": "MassDiskLensSmoothed",
            "type": MassDiskLensSmoothed,
            "params": { "density": 1.2, "radius": 3*AS, "endradius": 3.5*AS },
            "direct": False,
            "subpolynomial": True
        },
        {
            "name": "MultipleWendlandLens",
            "type": MultipleWendlandLens,
            "params": { 
                "phix": [
                    { "weight": 1, "scale": 100*AS, "position": [ 1*AS, 2*AS ] },
                    { "weight": 2, "scale": 15*AS, "position": [ -1*AS, -3*AS ] },
                ],
                "phiy": [
                    { "weight": 3, "scale": 30*AS, "position": [ 4*AS, 5*AS ] },
                    { "weight": 2, "scale": 45*AS, "position": [ -6*AS, 8*AS ] },
                    { "weight": 1, "scale": 65*AS, "position": [ -3*AS, 3*AS ] },
                ],
            },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "MultipleWendlandLens",
            "type": MultipleWendlandLens,
            "params": { 
                "points": [ 
                    [ -1*ANGLE_ARCSEC, 2*ANGLE_ARCSEC ],
                    [ 5*ANGLE_ARCSEC, 3*ANGLE_ARCSEC ],
                    [ 2*ANGLE_ARCSEC, -2*ANGLE_ARCSEC ],
                ],
                "angles": [ 
                    [ 2*ANGLE_ARCSEC, -2*ANGLE_ARCSEC ],
                    [ 5*ANGLE_ARCSEC, 3*ANGLE_ARCSEC ],
                    [ -1*ANGLE_ARCSEC, 2*ANGLE_ARCSEC ],
                ],
                "scale": 100*AS
            },
            "direct": False,
            "subpolynomial": False
        },
        { 
            "name": "NFWLens",
            "type": NFWLens,
            "params": { "rho_s": 1000,
                        "theta_s": 10*AS },
            "direct": True,
            "subpolynomial": False
        },
        { 
            "name": "EllipticNFWLens",
            "type": EllipticNFWLens,
            "params": { "rho_s": 1000,
                        "theta_s": 10*AS,
                        "q": 0.8 },
            "direct": True,
            "subpolynomial": False
        },
        { 
            "name": "SersicLens",
            "type": SersicLens,
            "params": { "index": 1.2,
                        "centraldensity": 2.1,
                        "scale": 5*AS },
            "direct": True,
            "subpolynomial": False
        },
        { 
            "name": "EllipticSersicLens",
            "type": EllipticSersicLens,
            "params": { "index": 1.2,
                        "centraldensity": 2.1,
                        "scale": 5*AS,
                        "q": 0.7 },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "ProfileLens",
            "type": ProfileLens,
            "params" : {
                "profile": [ 3, 2.5, 1.5, 0.1, 0],
                "radius": 10*AS
            },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "PIEMDLens",
            "type": PIEMDLens,
            "params": {
                "centraldensity": getCriticalDensity(Dd, 1.0, 0.8)*20,
                "coreradius": 0.2*ANGLE_ARCSEC,
                "scaleradius": 200*ANGLE_ARCSEC,
                "epsilon": 0.3
            },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "PIMDLens",
            "type": PIMDLens,
            "params": {
                "centraldensity": getCriticalDensity(Dd, 1.0, 0.8)*20,
                "coreradius": 0.2*ANGLE_ARCSEC,
                "scaleradius": 200*ANGLE_ARCSEC,
            },
            "direct": True,
            "subpolynomial": False
        },
        {
            "name": "HarmonicLens",
            "type": HarmonicLens,
            "params": {
                "sigma0": getCriticalDensity(Dd, 1.0, 0.8)*20,
                "k": 1/ANGLE_ARCSEC,
                "l": 2/ANGLE_ARCSEC,
                "phi_x": 1.2,
                "phi_y": 3.4
            },
            "direct": True,
            "subpolynomial": False
        }
]

for lensInfo in testLenses:
    t = lensInfo["type"]
    print("Checking lens {} ({})".format(lensInfo["name"], t))

    params = copy.deepcopy(lensInfo["params"])
    lens = t(Dd, params)

    if not lensInfo["subpolynomial"] and lensInfo["direct"]:
        obtainedParams = lens.getLensParameters()
        #pprint.pprint(obtainedParams)
        if obtainedParams == lensInfo["params"]:
            print("  Match")
        else:
            print("  NO MATCH")
    else:
        if lensInfo["subpolynomial"]:
            obtainedParams = lens.getPolynomialLensParameters()
            t = PolynomialMassProfileLens
        else:
            obtainedParams = lens.getLensParameters()

        newLens = t(Dd, obtainedParams)
        newParams = newLens.getLensParameters()

        if newParams == obtainedParams:
            print("  Match (indirect)")
        else:
            print("  NO MATCH (indirect)")
        
# DeflectionGridLens separately!
origParams = { 
    "bottomleft": [-10*AS, -10*AS], 
    "topright": [11*AS, 12*AS],
    "angles": np.array([
        [ [ -2*AS, -3*AS ], [ -1*AS, 2*AS] ],
        [ [ -3*AS, 5*AS ], [ 2*AS, -3*AS ]],
        [ [ 1*AS, 2*AS ], [ 3*AS, -2*AS ]],
    ], dtype=np.double)
}

lens = DeflectionGridLens(Dd, copy.deepcopy(origParams))
print("Checking lens {} ({})".format("DeflectionGridLens", type(lens)))
obtainedParams = lens.getLensParameters()
if (np.array_equal(obtainedParams["angles"], origParams["angles"]) and
   obtainedParams["bottomleft"] == origParams["bottomleft"] and
   obtainedParams["topright"] == origParams["topright"] and len(obtainedParams) == 3):
    print("  Match")
else:
    print("  NO MATCH")

# CompositeLens separately
origParams = [
    { "lens": PlummerLens(Dd, { "mass": 2e13*MASS_SUN, "width": 3*AS }),
      "factor": 0.8,
      "angle": 40,
      "x": -2*AS,
      "y": -4*AS
    },
    { "lens": GaussLens(Dd, { "mass": 1e13*MASS_SUN, "width": 4*AS }),
      "factor": 1.1,
      "angle": 3,
      "x": -4*AS,
      "y": 2*AS
    },
    { "lens": SquareLens(Dd, { "mass": 4e13*MASS_SUN, "width": 1*AS }),
      "factor": 0.3,
      "angle": 0.4,
      "x": 2*AS,
      "y": 1*AS
    },
]

lens = CompositeLens(Dd, origParams)
print("Checking lens {} ({})".format("CompositeLens", type(lens)))
obtainedParams = lens.getLensParameters()
isMatch = False
if len(obtainedParams) == len(origParams):
    for i in range(len(origParams)):
        if obtainedParams[i]["lens"].toBytes() == origParams[i]["lens"].toBytes():
            if obtainedParams[i]["factor"] == origParams[i]["factor"]:
                if obtainedParams[i]["angle"] == origParams[i]["angle"]:
                    if obtainedParams[i]["x"] == origParams[i]["x"]:
                        if obtainedParams[i]["y"] == origParams[i]["y"]:
                            if len(obtainedParams[i]) == 5:
                                isMatch = True
                            else:
                                print("Param component length error")
                        else:
                            print("Mismatch in y")
                    else:
                        print("Mismatch in x")
                else:
                    print("Mismatch in angle")
                    print(" Expected:", origParams[i]["angle"])
                    print(" Received:", obtainedParams[i]["angle"])
            else:
                print("Mismatch in factor")
                print(" Expected:", origParams[i]["factor"])
                print(" Received:", obtainedParams[i]["factor"])
        else:
            print("Mismatch in sublens bytes")
else:
    print("Param length error")

if isMatch:
    print("  Match")
else:
    print("  NO MATCH")
