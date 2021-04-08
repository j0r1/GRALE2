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

z_lens = 0.45
iws = inversion.InversionWorkSpace(z_lens, 150*ANGLE_ARCSEC, cosmology=Cosmology(0.7, 1.0, 0, 0))

iws.addImageDataToList(images.ImagesData.load("images1pointgroups.imgdata"), 2.5, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("null1.imgdata"), 2.5, "extendednullgrid")
iws.addImageDataToList(images.ImagesData.load("images2.imgdata"), 1.5, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("null2.imgdata"), 1.5, "extendednullgrid")

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

