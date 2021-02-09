# This is an optimization procedure, similar to the one in 
# "Lensing degeneracies and mass substructure"
# (https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.1772L/abstract)
# The other script also contains some information about the
# basic procedure.

from grale.all import *
import random

inversion.setDefaultInverter("mpi")

# These are the lens parameters that were used in making the model
cosm = cosmology.Cosmology(0.70, 0.27, 0, 0.73)
cosmology.setDefaultCosmology(cosm)

z_lens = 0.5
Dd = cosm.getAngularDiameterDistance(z_lens)

src1, src2 = images.ImagesData.load("src1.imgdat"), images.ImagesData.load("src2.imgdat")

# As explained in the other script, this is the MSD lens
# that works for src1, but not for the other one
baseLens = lenses.GravitationalLens.load("baselens.lensdata")

# The routine 'createMonopoleBasisFunctions' will call a callback
# for each grid cell. At the center of that cell, we'll monitor
# the density, to make sure it doesn't become negative. This ImagesData 
# instance will monitor those positions
densGrid = images.ImagesData(1, kappa=True)

def cellCenterCallback(ctr, cellSize, hasBasisFunction, densGrid):
    densGrid.addPoint(0, ctr, kappa=0)

# Here we create monopole basisfunctions based on a grid, with similar
# settings as in the article. The parameters also specifiy that we
# want to leave the images of src1 alone, while we do want overlap
# with the images of src2 (otherwise the deflection angles there
# cannot be affected)
#
# randomOffset is set to False to get a similar layout of basisfunctions
# as in the article; undoubtedly for better results you should set this
# to True and average several solutions.
basisFunctions = util.createMonopoleBasisFunctions([src1], Dd, 16, 80*ANGLE_ARCSEC*16/15, 
                        overlapNeeded=[src2], randomOffset=False,
                        cellCenterCallback=cellCenterCallback,
                        cellCenterCallbackState=densGrid
                        )

# This is not really the mass, but the rescaling routine in the GA
# needs something relevant. 
totalFakeMass = sum([ b["mass"] for b in basisFunctions ])

# For the inversion/optimization, we'll use two fitness measures:
# one to avoid negative densities, and one to try to restore the
# deflection angles at the images of source 2 to those of the original
# NSIE lens (these are stored in that file, see the other script)
iws = inversion.InversionWorkSpace(z_lens, 80*ANGLE_ARCSEC)
iws.addImageDataToList(densGrid, None, "kappathresholdmin")
iws.addImageDataToList(src2, None, "deflectionangles")
iws.addBasisFunctions(basisFunctions)

# Originally, in the 2012 article, a separate GA was used. The code
# now has become sufficiently general to optimize this with the
# default inversion algorithm
popSize = 512

lens, _, _ = iws.invertBasisFunctions(popSize,
    # Needed to make the GA look for a more appropriate scale factor
    massScale=totalFakeMass,
    fitnessObjectParameters={
        # We'll add an extra step in the convergence checks, that allows
        # some smaller modifications to the basis functions as well
        "convergence_smallmutation_sizes": np.array([-1.0,0.1,0.02]),
        "convergence_factors": np.array([0.1,0.05,0.01]),
        # This makes sure that the rescaling in the GA looks for the
        # best deflection angle fitness; otherwise the kappa threshold
        # (used to avoid negative densities) would be used.
        "scalepriority_deflectionangle": 0
    },
    baseLens=baseLens,
    allowNegativeValues=True,
)

# Note that these are just the corrections to the base lens, you'll need
# to add them using a CompositeLens for a full solution. If you run this
# script multiple times, you can first average the different monopoles.lensdata
# solutions, and then add that to the base lens.
lens.save("monopoles.lensdata")

