from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.plotutil as plotutil
from grale.cosmology import Cosmology
import grale.images as images
import numpy as np

# Write the RNG state, in case we want to reproduce the run exactly
# (note that the GRALE_DEBUG_SEED environment variable will need
# to be restored as well)
import random
print("RNG State:")
print(random.getstate())

V = lambda x, y: np.array([x,y], dtype=np.double)

renderers.setDefaultLensPlaneRenderer("mpi") # threads, mpi, opencl, None or a Renderer object
renderers.setDefaultMassRenderer("mpi") # threads, mpi, None, or a Renderer object
inversion.setDefaultInverter("mpi") # threads, mpi or an Inverter object
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

z_lens = 0.7
iws = inversion.InversionWorkSpace(z_lens, 60*ANGLE_ARCSEC, cosmology=Cosmology(0.7, 0.3, 0, 0.7))

iws.addImageDataToList(images.ImagesData.load("images_00.imgdata"), 2.74714, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_01.imgdata"), 3.07785, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_02.imgdata"), 2.85271, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_03.imgdata"), 3.32098, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_04.imgdata"), 3.1543 , "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_12.imgdata"), 3.00706, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_13.imgdata"), 3.28815, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_73.imgdata"), 3.3154 , "extendedimages")

# For a quick test to see if all code executes
# iws.setDefaultInversionArguments(maximumGenerations = 10)

iws.setUniformGrid(15)
lens1, fitness, fitdesc = iws.invert(512)
lens1.save("inv1.lensdata")

iws.setSubdivisionGrid(lens1, 300, 400)
lens2, fitness, fitdesc = iws.invert(512)
lens2.save("inv2.lensdata")

iws.setSubdivisionGrid(lens2, 500, 600)
lens3, fitness, fitdesc = iws.invert(512)
lens3.save("inv3.lensdata")

iws.setSubdivisionGrid(lens3, 700, 800)
lens4, fitness, fitdesc = iws.invert(512)
lens4.save("inv4.lensdata")

iws.setSubdivisionGrid(lens4, 900, 1000)
lens5, fitness, fitdesc = iws.invert(512)
lens5.save("inv5.lensdata")

