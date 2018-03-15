import grale.lenses as lenses
import grale.cosmology as cosmology
from grale.constants import *

# Create a cosmological model with dimensionless Hubble constant h=0.7,
# Omega_m = 0.3 and Omega_lambda = 0.7
cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)

# For a lens at z = 0.5, calculate the angular diameter distance
z_d = 0.5
D_d = cosm.getAngularDiameterDistance(z_d)

# Let's show what this angular diameter distance is, in units of Mpc
print("D_d = {:.2f} Mpc".format(D_d/DIST_MPC))

# Create a SIS lens with a velocity dispersion of 250 km/s
sis = lenses.SISLens(D_d, { "velocityDispersion": 250000 })

# Create a SIE lens with specific velocity dispersion and ellipticity
sie = lenses.SIELens(D_d, { "velocityDispersion": 200000, "ellipticity": 0.7 })

# Combine these models using a CompositeLens
combLens = lenses.CompositeLens(D_d, [
        { "lens": sis, "factor": 1,
          "x": -1*ANGLE_ARCSEC, "y": -2*ANGLE_ARCSEC, "angle": 0 },
        { "lens": sie, "factor": 1,
          "x": 2*ANGLE_ARCSEC, "y": 1.5*ANGLE_ARCSEC, "angle": 60 }
    ])


