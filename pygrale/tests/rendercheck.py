from grale.all import *
from grale.constants import *
from grale.renderers import *
import pprint
import numpy as np
import sys

class MyFeedback(feedback.Feedback):
    def onStatus(self, s):
        #print("STATUS", s)
        pass

    def onProgress(self, x):
        #sys.stdout.write("PROGRESS: {}    \r".format(x))
        pass

def check(name, lensModel, bottomLeft, topRight, numX, numY):
    print("=== Checking: {} ===".format(name))
    p = lensModel # used to be a plummer lens, hence the 'p'
    inputPoints = [ ]
    for y in range(numY):
        for x in range(numX):
            inputPoints.append([
                bottomLeft[0] + (topRight[0]-bottomLeft[0])/(numX-1) * x,
                bottomLeft[1] + (topRight[1]-bottomLeft[1])/(numY-1) * y
            ])

    thetas = np.array(inputPoints)
    angles = p.getAlphaVector(thetas)
    dens = p.getSurfaceMassDensity(thetas)

    f = MyFeedback()

    print("Grid based methods for deflections")
    for renderClass in [ ThreadsLensPlaneRenderer, OpenCLLensPlaneRenderer, MPILensPlaneRenderer ]:
        print("Using class:", renderClass)
        r = renderClass(feedbackObject = f)
        result = r.render(p.toBytes(), bottomLeft, topRight, numX, numY)
        result = np.frombuffer(result, "double").reshape([-1,5]) # for each point: ax, ay, axx, ayy, axy
        result = result[:,0:2]
        diff = angles/ANGLE_ARCSEC - result/ANGLE_ARCSEC
        print("Max diff:", diff.max(), "Min diff", diff.min())
        print("")
        
    print("Grid based methods for density")
    for renderClass in [ ThreadsMassDensityRenderer, MPIMassDensityRenderer ]:
        print("Using class:", renderClass)
        r = renderClass(feedbackObject = f)
        result = r.render(p.toBytes(), bottomLeft, topRight, numX, numY)
        result = np.frombuffer(result, "double")
        diff = dens/ANGLE_ARCSEC - result/ANGLE_ARCSEC
        print("Max diff:", diff.max(), "Min diff", diff.min())
        print("")

    print("Point vector based methods for deflections")
    for renderClass in [ ThreadsLensPlaneRenderer, OpenCLLensPlaneRenderer, MPILensPlaneRenderer ]:
        print("Using class:", renderClass)
        r = renderClass(feedbackObject = f)
        result = r.renderXYVector(p.toBytes(), thetas)
        result = np.frombuffer(result, "double").reshape([-1,5]) # for each point: ax, ay, axx, ayy, axy
        result = result[:,0:2]
        diff = angles/ANGLE_ARCSEC - result/ANGLE_ARCSEC
        print("Max diff:", diff.max(), "Min diff", diff.min())
        print("")

    print("Point vector based methods for density")
    for renderClass in [ ThreadsMassDensityRenderer, MPIMassDensityRenderer ]:
        print("Using class:", renderClass)
        r = renderClass(feedbackObject = f)
        result = r.renderXYVector(p.toBytes(), thetas)
        result = np.frombuffer(result, "double")
        diff = dens/ANGLE_ARCSEC - result/ANGLE_ARCSEC
        print("Max diff:", diff.max(), "Min diff", diff.min())
        print("")

    print("")

checkParameters = [
        {
            "name": "MultiplePlummerLens",
            "lensModel": lenses.MultiplePlummerLens(1000*DIST_MPC, [
                { 'mass': 2e12*MASS_SUN, 'width': 5*ANGLE_ARCSEC, 'x': -1*ANGLE_ARCSEC, 'y': 2*ANGLE_ARCSEC },
                { 'mass': 3e12*MASS_SUN, 'width': 4*ANGLE_ARCSEC, 'x': 2*ANGLE_ARCSEC, 'y': -1*ANGLE_ARCSEC },
                ]),
            "bottomLeft": (-9*ANGLE_ARCSEC, -10*ANGLE_ARCSEC ),
            "topRight": (12*ANGLE_ARCSEC, 8*ANGLE_ARCSEC ),
            "numX": 131,
            "numY": 173,
        },
        {
            "name": "PlummerLens",
            "lensModel": lenses.PlummerLens(1000*DIST_MPC, { 'mass': 5e12*MASS_SUN, 'width': 5*ANGLE_ARCSEC }),
            "bottomLeft": (-9*ANGLE_ARCSEC, -10*ANGLE_ARCSEC ),
            "topRight": (12*ANGLE_ARCSEC, 8*ANGLE_ARCSEC ),
            "numX": 131,
            "numY": 173,
        },
        {
            "name": "SISLens",
            "lensModel": lenses.SISLens(1000*DIST_MPC, { 'velocityDispersion': 400000 }),
            "bottomLeft": (-9*ANGLE_ARCSEC, -10*ANGLE_ARCSEC ),
            "topRight": (12*ANGLE_ARCSEC, 8*ANGLE_ARCSEC ),
            "numX": 131,
            "numY": 173,
        },
        {
            "name": "SquareLens",
            "lensModel": lenses.SquareLens(1000*DIST_MPC, { 'mass': 5e12*MASS_SUN, 'width':2*ANGLE_ARCSEC }),
            "bottomLeft": (-9*ANGLE_ARCSEC, -10*ANGLE_ARCSEC ),
            "topRight": (12*ANGLE_ARCSEC, 8*ANGLE_ARCSEC ),
            "numX": 131,
            "numY": 173,
        },
        {
            "name": "MassSheetLens",
            "lensModel": lenses.MassSheetLens(1000*DIST_MPC, { 'Dds': 0.8, 'Ds': 1.0 }),
            "bottomLeft": (-9*ANGLE_ARCSEC, -10*ANGLE_ARCSEC ),
            "topRight": (12*ANGLE_ARCSEC, 8*ANGLE_ARCSEC ),
            "numX": 131,
            "numY": 173,
        },
        {
            "name": "CompositeLens",
            "lensModel": lenses.CompositeLens(1000*DIST_MPC, [
                {
                    "lens": lenses.CompositeLens(1000*DIST_MPC, [
                        { "lens": lenses.PlummerLens(1000*DIST_MPC, {"mass": 5e12*MASS_SUN, "width": 4*ANGLE_ARCSEC}), "factor": 1.01, "x": 2.5*ANGLE_ARCSEC, "y": -2.1*ANGLE_ARCSEC, "angle": 0 },
                        { "lens": lenses.SISLens(1000*DIST_MPC, {"velocityDispersion": 300000 }), "factor": 1.1, "x": -2*ANGLE_ARCSEC, "y": 1*ANGLE_ARCSEC, "angle": 0 },

                    ]), "x": 6*ANGLE_ARCSEC, "y": 5*ANGLE_ARCSEC, "factor": 0.9, "angle": 60
                },
                {
                    "lens": lenses.CompositeLens(1000*DIST_MPC, [
                        { "lens": lenses.SquareLens(1000*DIST_MPC, {"mass": 2e12*MASS_SUN, "width": 3*ANGLE_ARCSEC}), "factor": .3, "x": 2.5*ANGLE_ARCSEC, "y": -2.1*ANGLE_ARCSEC, "angle": 0 },
                        {
                            "lens": lenses.MultiplePlummerLens(1000*DIST_MPC, [
                                { "mass": 1e12*MASS_SUN, "width": 2*ANGLE_ARCSEC, "x": 1*ANGLE_ARCSEC, "y": 0.5*ANGLE_ARCSEC },
                                { "mass": 1e12*MASS_SUN, "width": 1.5*ANGLE_ARCSEC, "x": -2*ANGLE_ARCSEC, "y": -0.75*ANGLE_ARCSEC },
                            ]), "factor": 0.4, "x": -2*ANGLE_ARCSEC, "y": 1*ANGLE_ARCSEC, "angle": 30 },
                    ]), "x": -4*ANGLE_ARCSEC, "y": -3*ANGLE_ARCSEC, "factor": 1.3, "angle": -40
                },
                { "lens": lenses.MassSheetLens(1000*DIST_MPC, { 'Dds': 0.8, 'Ds': 1.0 }), "factor": 0.1, "x": 1*ANGLE_ARCSEC, "y": 2*ANGLE_ARCSEC, "angle": 10 },
            ]),
            "bottomLeft": (-11*ANGLE_ARCSEC, -10*ANGLE_ARCSEC ),
            "topRight": (12*ANGLE_ARCSEC, 11*ANGLE_ARCSEC ),
            "numX": 131,
            "numY": 173,
        }
]

def main():
    if len(sys.argv) == 8:
        check(
                "Command line model",
                lenses.GravitationalLens.load(sys.argv[1]),
                V(float(sys.argv[2]), float(sys.argv[3]))*ANGLE_ARCSEC,
                V(float(sys.argv[4]), float(sys.argv[5]))*ANGLE_ARCSEC,
                int(sys.argv[6]),
                int(sys.argv[7]) )
    elif len(sys.argv) == 1:
        for params in checkParameters:
            check(**params)
    else:
        print("Either specify no parameters to run default checks, or specify lensmodel, bottomleft, topright, numx and numy")

if __name__ == "__main__":
    main()
