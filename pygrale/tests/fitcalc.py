import grale.inversion as inversion
import grale.lenses as lenses
import grale.feedback as feedback
import grale.cosmology as cosmology
import grale.plotutil as plotutil
import grale.images as images
from grale.constants import *
import numpy as np

feedback.setDefaultFeedback("none")

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)
D = cosm.getAngularDiameterDistance

def main():
    zd = 0.5
    lens = lenses.PlummerLens(D(zd), { "mass": 1e14*MASS_SUN, "width": 5*ANGLE_ARCSEC })
    li = plotutil.LensInfo(lens, size=15*ANGLE_ARCSEC, zd=zd)

    iws = inversion.InversionWorkSpace(zd, 20*ANGLE_ARCSEC)
    
    num = 5
    print(f"Creating {num} sources")
    while num > 0:
        zs = np.random.uniform(zd, 6)
        beta = np.random.uniform(-1*ANGLE_ARCSEC, 1*ANGLE_ARCSEC, size=(2,))
        
        li.setSourceRedshift(zs)
        ip = li.getImagePlane()
        thetas = ip.traceBeta(beta)
        if len(thetas) > 1:
            num -= 1
            print(f"Found image, {num} to go")

            img = images.ImagesData(len(thetas))
            for i in range(len(thetas)):
                img.addPoint(i, thetas[i])

            iws.addImageDataToList(img, zs, "pointimages")

    print("Fitness info")
    print("------------")
    print(iws.calculateFitness(None))
    print("")


    imperfectLens = lenses.CompositeLens(D(zd), [
        { "lens": lens, "factor": 0.999, "x": 0, "y": 0, "angle": 0 } ])

    for l, info in [ (lens, "True lens"), (imperfectLens, "Imperfect lens") ]:
        s = f"Using: {info}"
        print(s)
        print("-" * len(s))
        print("Using lens itself:")
        print(iws.calculateFitness(l))
        print("")

        print("Using backprojected images:")
        bpImages = iws.backProject(l)
        print(iws.calculateFitness(bpImages))
        print("")
        

if __name__ == "__main__":
    main()
