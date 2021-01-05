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

strongSize = 200*ANGLE_ARCSEC
weakSize = 30*ANGLE_ARCMIN 
weakSubDiv = 48
weakThreshold = 0.1 # Threshold for |1-kappa|
sheetType = "nosheet" # or "genome"

iws = inversion.InversionWorkSpace(z_lens, strongSize)

# Add the SL data
# Note: not using null space
for i in images.readInputImagesFile("images.txt", True):
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")

iws.addImageDataToList(images.ImagesData.load("ellipt_48x48_exact_z1.imgdata"), 1, "sheardata", { "threshold": weakThreshold })
iws.addImageDataToList(images.ImagesData.load("ellipt_48x48_exact_z2.imgdata"), 2, "sheardata", { "threshold": weakThreshold })
iws.addImageDataToList(images.ImagesData.load("ellipt_48x48_exact_z4.imgdata"), 4, "sheardata", { "threshold": weakThreshold })

iws.setDefaultInversionArguments(sheetSearch = sheetType) # add maximumGenerations = 2 to test if script works

# This is something that can be used below:
# The genomes will only describe the mass in each basisfunction up to a scale
# factor. The mass in the strong lensing region can be estimated roughly from
# the image positions, and the usual strategy is to look for an appropriate 
# scale factor within a certain range around this estimate.
# If we're also including basis functions in the WL area, by default their
# mass is also taken account to estimate the total lens mass, which can
# differ a lot from the SL mass estimate. There are two strategies to handle
# this:
#  1. increase the range in which the scale factor is sought
#  2. don't count the WL basis function masses in the rescaling operation, so
#     that that's only based on the SL mass
#
# The first can be done by setting 'massScaleSearchType' to 'noshearch' in an invert
# function. The second is what the function below is for: We'll just call
# the default implementation, which returns a basis function, as well as the
# mass that is counted for the rescaling operation. By setting this to something
# very small, the basis function's mass is ignored when determining an 
# appropriate scale factor.
# 
# Note that there is now a convenience function setStrongAndWeakBasisFunctions
# in InversionWorkSpace that does this internally.

def myLensModelFunction(operation, operationInfo, parameters):
    r = inversion.defaultLensModelFunction(operation, operationInfo, parameters)
    if operation == "add":
        return (r[0], 0.1) # Set the mass that's counted in determining the scale factor in the GA to (almost) zero
    return r

def setBasisFunctions():
    iws.clearBasisFunctions()
    iws.addBasisFunctionsBasedOnCurrentGrid() # For the SL data

    # Create a new, uniform grid for the WL part
    iws.setUniformGrid(weakSubDiv, regionSize = weakSize)
    # Here we use our own function instead of the default (see above). The basis
    # functions are initialized to yield a total mass of 1e15 M_sun
    iws.addBasisFunctionsBasedOnCurrentGrid(myLensModelFunction, { "totalmass": 1e15*MASS_SUN })

iws.setUniformGrid(15)
setBasisFunctions()
lens1, fitness1, fitdesc1 = iws.invertBasisFunctions(512)
lens1.save("inv1.lensdata")

iws.setSubdivisionGrid(lens1, 300, 400)
setBasisFunctions()
lens2, fitness2, fitdesc2 = iws.invertBasisFunctions(512)
lens2.save("inv2.lensdata")

iws.setSubdivisionGrid(lens2, 500, 600)
setBasisFunctions()
lens3, fitness3, fitdesc3 = iws.invertBasisFunctions(512)
lens3.save("inv3.lensdata")

iws.setSubdivisionGrid(lens3, 700, 800)
setBasisFunctions()
lens4, fitness4, fitdesc4 = iws.invertBasisFunctions(512)
lens4.save("inv4.lensdata")

iws.setSubdivisionGrid(lens4, 900, 1000)
setBasisFunctions()
lens5, fitness5, fitdesc5 = iws.invertBasisFunctions(512)
lens5.save("inv5.lensdata")

