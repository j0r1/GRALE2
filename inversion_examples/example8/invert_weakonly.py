from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.images as images
import grale.lenses as lenses
import grale.cosmology as cosmology

# Write the RNG state, in case we want to reproduce the run exactly
# (note that the GRALE_DEBUG_SEED environment variable will need
# to be restored as well)
import random
print("RNG State:")
print(random.getstate())

renderers.setDefaultLensPlaneRenderer("mpi") # threads, mpi, opencl, None or a Renderer object
renderers.setDefaultMassRenderer("mpi") # threads, mpi, None, or a Renderer object
inversion.setDefaultInverter("mpi") # singlecore, mpi, localcs, mpics, or an Inverter object

z_lens = 0.4
cosm = cosmology.Cosmology(0.7, 0.27, 0, 0.73)
cosmology.setDefaultCosmology(cosm)

weakSize = 30*ANGLE_ARCMIN 
weakSubDiv = 48
weakThreshold = 0.1 # Threshold for |1-kappa|
sheetType = "nosheet" # or "genome"

iws = inversion.InversionWorkSpace(z_lens, weakSize) # Now using weakSize for the default grid size
iws.addImageDataToList(images.ImagesData.load("ellipt_48x48_exact_z1.imgdata"), 1, "sheardata", { "threshold": weakThreshold })
iws.addImageDataToList(images.ImagesData.load("ellipt_48x48_exact_z2.imgdata"), 2, "sheardata", { "threshold": weakThreshold })
iws.addImageDataToList(images.ImagesData.load("ellipt_48x48_exact_z4.imgdata"), 4, "sheardata", { "threshold": weakThreshold })

iws.setDefaultInversionArguments(sheetSearch = sheetType) # add maximumGenerations = 2 to test if script works

iws.setUniformGrid(weakSubDiv)

# We don't have SL info to estimate a mass scale, so we'll specify one. This
# will be used to set the masses of the basis functions for each grid cell.

# For WL only inversions, the extra search for a mass scale does not really
# provide much benefit, just overhead. For this reason, it is set to
# "nosearch" here. 

lens, fitness, fitdesc = iws.invert(512, massScale=5e15*MASS_SUN, massScaleSearchType="nosearch")
lens.save("inv_weakonly.lensdata")
