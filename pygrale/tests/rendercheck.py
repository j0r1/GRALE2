from __future__ import print_function
import grale.lenses as lenses
import grale.images as images
import grale.feedback as feedback
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

def main():
    p = lenses.PlummerLens(1000*DIST_MPC, { 'mass': 5e12*MASS_SUN, 'width': 5*ANGLE_ARCSEC }) 

    bottomLeft = (-9*ANGLE_ARCSEC, -10*ANGLE_ARCSEC )
    topRight = (12*ANGLE_ARCSEC, 8*ANGLE_ARCSEC )
    numX = 131
    numY = 173

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

if __name__ == "__main__":
    main()
