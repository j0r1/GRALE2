import grale.lenses as lenses
import grale.multiplane as multiplane
import grale.cosmology as cosmology
from grale.constants import *
import numpy as np

V = lambda x, y: np.array([x, y], dtype=np.double)

h = np.random.uniform(0.6, 0.8)
Wm = np.random.uniform(0.25, 0.35)
Wr = np.random.uniform(0, 0.1)
if np.random.random() < 0.5:
    Wv = 1-Wr-Wm
else:
    Wv = np.random.uniform(0.6, 0.8)
w = np.random.uniform(-1.1, -0.9)

cosm = cosmology.Cosmology(h, Wm, Wr, Wv, w)
D = cosm.getAngularDiameterDistance
print(f"const Cosmology cosm({h}, {Wm}, {Wr}, {Wv}, {w});")

cosmology.setDefaultCosmology(cosm)
numLensPlanes = np.random.randint(2, 5)

lensRedShifts = np.random.uniform(0.4, 2.5, (numLensPlanes,))
# TODO: This does not seem to work if not sorted!
lensRedShifts = sorted(lensRedShifts)
lensplanes = []

def getNextPlummerParams():
    mass = np.random.uniform(1e12*MASS_SUN, 1e13*MASS_SUN)
    width = np.random.uniform(2, 5)*ANGLE_ARCSEC
    x = np.random.uniform(-10, 10)*ANGLE_ARCSEC
    y = np.random.uniform(-10, 10)*ANGLE_ARCSEC
    return mass, width, x, y

# dbgPlummers = [ (1e14*MASS_SUN, 5*ANGLE_ARCSEC, 0, 0),
#                 (5e13*MASS_SUN, 3*ANGLE_ARCSEC, 1*ANGLE_ARCSEC, -1*ANGLE_ARCSEC),
# ]
# def getNextPlummerParams():
#     return dbgPlummers.pop(0)

def getNextTheta():
    return V(np.random.uniform(-40, 40), np.random.uniform(-40, 40))*ANGLE_ARCSEC

# dbgThetas = [ V(10,8)*ANGLE_ARCSEC, V(-12,2)*ANGLE_ARCSEC ]
# def getNextTheta():
#     return dbgThetas.pop(0)

print("const vector<pair<double, vector<PlummerLensInfo>>> redshiftsAndPlummers = {");
totalMass = 0
for zd in lensRedShifts:
    print(f"  {{ {zd}, {{ ")
    Dd = D(zd)
    numPlummers = np.random.randint(5, 10)

    params = []
    for i in range(numPlummers):
        mass, width, x, y = getNextPlummerParams()
        totalMass += mass

        props = {
            "lens": lenses.PlummerLens(Dd, {
                "mass": mass,
                "width": width
            }),
            "x": x,
            "y": y,
            "angle": 0,
            "factor": 1
        }
        params.append(props)
        print(f"    {{ {mass}, {width}, {{ {x}, {y} }} }},")

    print("  } },")
    # Just store the parameters at this point, we'll set the factor below
    # to get an approximate target mass
    lensplanes.append({ "zd": zd, "Dd": Dd, "params": params})

print("};")

imgSep = 30*ANGLE_ARCSEC
DdMin = min([p["Dd"] for p in lensplanes])
massForImgSep = imgSep**2 * DdMin * SPEED_C**2/(4*CONST_G)
print(f"// massForImgSep: {massForImgSep/MASS_SUN:g}")
print(f"// currentTotalMass: {totalMass/MASS_SUN:g}")
extraFactor = massForImgSep*np.random.uniform(0.95, 1.05)/totalMass

print(f"const double scaleFactor = {extraFactor}; // This is the scale factor for which we should find a match")

# Adjust individual plummer factors
for x in lensplanes:
    params = x["params"]
    for p in params:
        p["factor"] = extraFactor

lensesAndRedshifts = [ (lenses.CompositeLens(x["Dd"], x["params"]), x["zd"]) for x in lensplanes ]
# The grid dimensions are not used for this test
lensPlane = multiplane.MultiLensPlane(lensesAndRedshifts, V(-50,-50)*ANGLE_ARCSEC, V(50, 50)*ANGLE_ARCSEC, 64, 64)

numSources = np.random.randint(10, 20)
sourceRedshifts = np.random.uniform(min(lensRedShifts)+0.1, 5.0, numSources)
# TODO: this doesn't seem to work if not sorted!
# sourceRedshifts = sorted(sourceRedshifts)

print("const vector<pair<double,vector<pair<Vector2Dd, Vector2Dd>>>> thetaBetaMappings {")
for zs in sourceRedshifts:
    imgPlane = multiplane.MultiImagePlane(lensPlane, zs)

    numPoints = np.random.randint(3, 7)

    thetas = np.array([ getNextTheta() for i in range(numPoints) ])
    betas = imgPlane.traceTheta(thetas)

    print(f"  {{ {zs}, {{")
    for i in range(len(thetas)):
        print(f"    {{ {{ {thetas[i][0]}, {thetas[i][1]}  }},  {{ {betas[i][0]}, {betas[i][1]}  }} }},")
    print("  } }, ")
print("};")
