from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.images as images
import grale.lenses as lenses
import grale.cosmology as cosmology

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
# We don't have SL info to estimate a mass scale, so we'll specify one. If you're
# more uncertain about this, you can also specify wideSearch=True
lens, fitness, fitdesc = iws.invert(512, massScale=5e15*MASS_SUN)
lens.save("inv_weakonly.lensdata")
