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


# Set the source redshift to 4
zs = 4.0
lensInfo = LI( [(movedSis, z1), (movedSie, z2)], size=6*ANGLE_ARCSEC, zs=zs)

# Create a circular source at a position
V = lambda x,y: np.array([x,y], dtype=np.double)
src = images.CircularSource(V(0.35,0.1)*ANGLE_ARCSEC, 0.05*ANGLE_ARCSEC, 1, True)

# And plot the situation
plotutil.plotImagePlane(lensInfo, [src])
plt.show()
