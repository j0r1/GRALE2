import grale.images as images
import grale.lenses as lenses
import grale.cosmology as cosmology
import grale.plotutil as plotutil
from grale.constants import *
import numpy as np
import matplotlib.pyplot as plt

# Set a default cosmological model and create an alias for the
# angular diameter function
cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)
Z = cosm.getAngularDiameterDistance

# Use a default plot unit of one arcsec
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

# A shorter name for the LensInfo class
LI = plotutil.LensInfo

# Create a SIS lens with a velocity dispersion of 250 km/s and
# set it at a specific position
z1 = 0.5
sis = lenses.SISLens(Z(z1), { "velocityDispersion": 250000 })
movedSis = lenses.CompositeLens(Z(z1), [
        { "lens": sis, "factor": 1,
          "x": -0.25*ANGLE_ARCSEC, "y": -0.5*ANGLE_ARCSEC, "angle": 0 }
])

# Create a SIE lens with specific velocity dispersion and ellipticity
# at different redshift
z2 = 1.0
sie = lenses.SIELens(Z(z2), { "velocityDispersion": 200000, "ellipticity": 0.7 })
movedSie = lenses.CompositeLens(Z(z2), [
        { "lens": sie, "factor": 1,
          "x": 1*ANGLE_ARCSEC, "y": 0.75*ANGLE_ARCSEC, "angle": 60 }
    ])

# A shorter name for the DensInfo class
DI = plotutil.DensInfo

# Create a LensInfo instance for the SIS lens, which we'll use to
# calculate the lens effect on the shape of the SIE
# Here, we specify redshifts instead of angular diameter distances.
lensInfoSIS = LI(movedSis, zd=z1, zs=z2, size=6*ANGLE_ARCSEC)
# For the SIE, we'll only use this to get the shape of the mass
# distribution (as pixels)
lensInfoSIE = LI(movedSie, size=6*ANGLE_ARCSEC)

# Get the densities for both lenses
sisDens = lensInfoSIS.getDensityPixels()
sieDens = lensInfoSIE.getDensityPixels()

# This describes the combination of the densities, without taking the lens
# effect by the foreground SIS into account
comboUnlensed = DI(sieDens + sisDens, **lensInfoSIS.getArea(), isPixels = True)

# To calculate the distortion on the SIE, we'll use its density as a source 
# shape. Because the DiscreteSource expects data in the order of an image (top left
# pixel first), we need to flip it.
sieSrc = images.DiscreteSource(np.flipud(sieDens), 6*ANGLE_ARCSEC, 6*ANGLE_ARCSEC, [0,0])

plt.figure(figsize=(10,10))
plt.subplot(2,2,1)

# First, plot the SIE density itself
plotutil.plotDensity(lensInfoSIE, vmax=10)
plt.gca().set_title("Unlensed SIE at z={}".format(z2))

# We're going to store the lensed image by using this callback as the
# 'processRenderPixels' value. We'll also rescale the brightness a bit, since the
# plotImagePlane routine expects values between 0 and 1. We'll also clip the final
# result, for a cleaner result.
lensedSIE = [ None ]
def f(x):
    if lensedSIE[0] is None:
        # 'x' contains RGB value, which are 3 same values for a gray scale image
        # We're only going to store one of them, to combine it later with the
        # SIS density (which is also a single value per pixel)
        lensedSIE[0] = x[:,:,0].reshape(x.shape[0:2])
    return np.clip(x/10, 0, 1)

plt.subplot(2,2,2)

# First we call the routine without the sources, so we'll only store the lensed
# SIE shape
plotutil.plotImagePlane(lensInfoSIS, [sieSrc], plotSources=False, processRenderPixels = f)

# In the second call, the pixels will no longer be stored in our function
# but the rescaling will still be done.
plotutil.plotImagePlane(lensInfoSIS, [sieSrc], processRenderPixels = f)
plt.gca().set_title("Lensed SIE by SIS at z={}".format(z1))

plt.subplot(2,2,3)
plotutil.plotDensity(comboUnlensed, vmax = 10)
plt.gca().set_title("Both, without lensing")

# Combine the lensed SIE shape with the SIS shape
comboLensed = DI(lensedSIE[0] + sisDens, **lensInfoSIS.getArea(), isPixels = True)
plt.subplot(2,2,4)
plotutil.plotDensity(comboLensed, vmax = 10)
plt.gca().set_title("Both, with lensing")
plt.show()

