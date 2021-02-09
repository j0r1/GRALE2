# Here, we create some input for the inversion/optimization that
# is similar to the one from "Lensing degeneracies and mass substructure"
# (https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.1772L/abstract)
#
# Basically, we create an NSIE lens and lens two sources. A simple
# mass-sheet degenerate (MSD) lens is created for the redshift of
# the first source, which then doesn't work for the second. In the
# optimization script, we'll try to restore the deflection angles
# at the images of the second source to those of the original NSIE
# lens, while leaving the deflection angles at the images of the
# first source alone. If that is successful, the deflection angles
# for source 1 are still that of the MSD version, while the deflection
# angles for source 2 are those of the original NSIE lens. So the
# images of the first source will correspond to a rescaled source,
# while those of the second source won't.

from grale.all import *

cosm = cosmology.Cosmology(0.7, 0.27, 0.0, 0.73)
cosmology.setDefaultCosmology(cosm)
D = cosm.getAngularDiameterDistance

zd = 0.5
zs1 = 1.2
zs2 = 1.8

src1Params = {
    "position": V(10, 5)*ANGLE_ARCSEC,
    "halfAxis": 1.1*ANGLE_ARCSEC,
    "eccentricity": 0.6,
    "angle": 30,
}

src2Params = {
    "position": V(2, -4.5)*ANGLE_ARCSEC,
    "halfAxis": 0.8*ANGLE_ARCSEC,
    "eccentricity": 0.4,
    "angle": 110,
}

l = lenses.NSIELens(D(zd), { "velocityDispersion": 1300000, "ellipticity": 0.8, "coreRadius": 1.5*ANGLE_ARCSEC })
li = plotutil.LensInfo(l, size=80*ANGLE_ARCSEC, zd=zd, numxy=511)

li.setSourceRedshift(zs1)
ip = li.getImagePlane()
src1ImagePixels = ip.segment(ip.renderImages([images.EllipticalSource(**src1Params)]))

li.setSourceRedshift(zs2)
ip = li.getImagePlane()
src2ImagePixels = ip.segment(ip.renderImages([images.EllipticalSource(**src2Params)]))

# We also store the deflection angles at the image points. We're not actually
# using those of the first source, but their presence doesn't cause a problem
# either
def pixelsToImagesData(pix):
    imgDat = images.ImagesData(len(pix), alpha=True)
    for i, img in zip(range(len(pix)), pix):
        for p in img:
            imgDat.addPoint(i, p, alpha=l.getAlphaVector(p)/ANGLE_ARCSEC)
    return imgDat

imgDat1 = pixelsToImagesData(src1ImagePixels)
imgDat2 = pixelsToImagesData(src2ImagePixels)

# These are the images that are used in the optimization procedure
imgDat1.save("src1.imgdat")
imgDat2.save("src2.imgdat")

lambd = 0.75

# Construct mass sheet degeneracy for source1
sheet = lenses.MassSheetLens(D(zd), { "Ds": D(zs1), "Dds": D(zd,zs1)})

degen = lenses.CompositeLens(D(zd), [
    { "lens": l, "factor": lambd, "x": 0, "y": 0, "angle": 0 },
    { "lens": sheet, "factor": (1-lambd), "x": 0, "y": 0, "angle": 0}
])

# This MSD lens works for imgDat1, but not for imgDat2. The optimization
# script will try to correct that.
degen.save("baselens.lensdata")

