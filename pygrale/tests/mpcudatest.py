import grale.multiplanecuda as multiplanecuda
import grale.multiplane as multiplane
import grale.cosmology as cosmology
import grale.lenses as lenses
import grale.plotutil as plotutil
import grale.renderers as renderers
import grale.images as images
import grale.feedback as feedback
from grale.constants import *
import matplotlib.pyplot as plt
import numpy as np
import pprint

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)
D = cosm.getAngularDiameterDistance

V = lambda x,y : np.array([x, y], dtype=np.double)
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)
feedback.setDefaultFeedback("none")

zd1, zd2 = 0.5, 1.0 
lens1 = lenses.CompositeLens(D(zd1), [
    { "factor": 0.7, "x": 1.0*ANGLE_ARCSEC, "y": -0.5*ANGLE_ARCSEC, "angle": 0 ,
      "lens": lenses.PlummerLens(D(zd1), { "mass": 1e14*MASS_SUN, "width": 4*ANGLE_ARCSEC } )},
    ])

lens2 = lenses.MultiplePlummerLens(D(zd2), [
    { "mass": 1e14*MASS_SUN, "width": 3*ANGLE_ARCSEC, "x": 0, "y": 2*ANGLE_ARCSEC },
    ])

lensesAndRedshifts = [
        (lens1, zd1),
        (lens2, zd2)
        ]

HW = 40*ANGLE_ARCSEC
li = plotutil.LensInfo(lensesAndRedshifts, size=2*HW, numxy=1023)
src1 = images.CircularSource(V(0.4,0.6)*ANGLE_ARCSEC, 0.5*ANGLE_ARCSEC)
zs1 = 0.9

src2 = images.CircularSource(V(10,10)*ANGLE_ARCSEC, 1.5*ANGLE_ARCSEC)
zs2 = 1.9
plotutil.plotImagePlane(li, [{"shape": src1, "z": zs1}, {"shape": src2, "z": zs2}], plotCaustics=False, plotCriticalLines=False)
plt.show()

thetasAndSourceRedshifts = [ ]
for zs, src in [ (zs1, src1), (zs2, src2)]:
    li.setSourceRedshift(zs)
    ip = li.getImagePlane()
    segs = ip.segment(ip.renderImages([src1]))
    thetas = np.concatenate(tuple(segs))
    thetasAndSourceRedshifts.append((thetas, zs))

mpCuda = multiplanecuda.MultiPlaneCUDA(lensesAndRedshifts, thetasAndSourceRedshifts)
print("Plummer parameters")
pprint.pprint(mpCuda.getPlummerParameters())

for factors in [ [[0.9],[1.1]], [[1],[1]] ]:
    mpCuda.calculateSourcePositions(factors)
    print("Showing results for", factors)
    print("(if not all 1, then the results from the CUDA version and the CPU version will differ,")
    print(" the CPU version does not take any factors into account)")
    print("")

    for i, thetas, zs, in zip(range(2), *zip(*thetasAndSourceRedshifts)):
        srcPos = mpCuda.getSourcePositions(i)/ANGLE_ARCSEC
        plt.plot(srcPos[:,0], srcPos[:,1], '.', label='CUDA')

        li.setSourceRedshift(zs)
        ip = li.getImagePlane()
        srcPos = ip.traceTheta(thetas)/ANGLE_ARCSEC
        plt.plot(srcPos[:,0], srcPos[:,1], '+', label='CPU')
    
        plt.gca().set_title(f"factors = {factors}, source idx = {i}")
        plt.gca().legend()
        plt.show()

