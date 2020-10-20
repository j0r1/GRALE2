from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.plotutil as plotutil
from grale.cosmology import Cosmology
import grale.images as images
import grale.lenses as lenses
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
inversion.setDefaultInverter("mpi") # singlecore, mpi, localcs, mpics, or an Inverter object
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

z_lens = 0.7
iws = inversion.InversionWorkSpace(z_lens, 60*ANGLE_ARCSEC, cosmology=Cosmology(0.7, 0.3, 0, 0.7))

iws.addImageDataToList(images.ImagesData.load("images_00.imgdata"), 2.74714, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_01.imgdata"), 3.07785, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_02.imgdata"), 2.85271, "extendedimages")

iws.addImageDataToList(images.ImagesData.load("images_03_td.imgdata"), 3.32098, "extendedimages")
#iws.addImageDataToList(images.ImagesData.load("images_03_td.imgdata"), 3.32098, "extendedimages", { "timedelay": False })

iws.addImageDataToList(images.ImagesData.load("images_04.imgdata"), 3.1543 , "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_12.imgdata"), 3.00706, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_13.imgdata"), 3.28815, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("images_73.imgdata"), 3.3154 , "extendedimages")

#iws.setDefaultInversionArguments(maximumGenerations = 10)
#iws.setDefaultInversionArguments(fitnessObjectParameters = { 'priority_timedelay':0 }, maximumGenerations = 10)
iws.setDefaultInversionArguments(geneticAlgorithmParameters = { "selectionpressure": 2 }, sheetSearch = "genome")
#iws.setDefaultInversionArguments(geneticAlgorithmParameters = { "selectionpressure": 2 }, sheetSearch = "genome", maximumGenerations=10)

iws.setUniformGrid(15)
lens1, fitness1, fitdesc = iws.invert(512)
lens1.save("inv1.lensdata")
baseLens, baseFitness = lens1, fitness1

iws.setSubdivisionGrid(lens1, 300, 400)
lens2, fitness2, fitdesc = iws.invert(512)
lens2.save("inv2.lensdata")
if fitness2 < baseFitness:
    baseLens, baseFitness = lens2, fitness2

iws.setSubdivisionGrid(lens2, 500, 600)
lens3, fitness3, fitdesc = iws.invert(512)
lens3.save("inv3.lensdata")
if fitness3 < baseFitness:
    baseLens, baseFitness = lens3, fitness3

iws.setSubdivisionGrid(lens3, 700, 800)
lens4, fitness4, fitdesc = iws.invert(512)
lens4.save("inv4.lensdata")
if fitness4 < baseFitness:
    baseLens, baseFitness = lens4, fitness4

iws.setSubdivisionGrid(lens4, 900, 1000)
lens5, fitness5, fitdesc = iws.invert(512)
lens5.save("inv5.lensdata")
if fitness5 < baseFitness:
    baseLens, baseFitness = lens5, fitness5

lensInfo = plotutil.LensInfo(baseLens, size=60*ANGLE_ARCSEC)
corrMassScale = lensInfo.getIntegratedMass()/10.0
print("Mass scale for corrections: {:g} solar masses".format(corrMassScale/MASS_SUN))

iws.setUniformGrid(64)
corrections, corrFitness, fitdesc = iws.invert(512, baseLens = baseLens, allowNegativeValues = True, massScale = corrMassScale)
corrections.save("corrections.lensdata")

correctedLens = lenses.CompositeLens(baseLens.getLensDistance(), [
                                     { "factor": 1, "x": 0, "y": 0, "angle": 0, "lens": baseLens },
                                     { "factor": 1, "x": 0, "y": 0, "angle": 0, "lens": corrections } ])
correctedLens.save("correctedlens.lensdata")

print(fitness1)
print(fitness2)
print(fitness3)
print(fitness4)
print(fitness5)
print(baseFitness)
print(corrFitness)

def compareTimeDelays(l):
    zd = 0.7
    zs = 3.32098
    imgTd = images.ImagesData.load("images_03_td.imgdata");
    cosm = iws.getCosmology()
    Ds = cosm.getAngularDiameterDistance(zs)
    Dds = cosm.getAngularDiameterDistance(zd, zs)

    betaAvg = V(0,0)
    thetas = [ ]
    print("Time delays stored in the input file:")
    for i in range(imgTd.getNumberOfTimeDelays()):
        imgNum, ptNum, timeDelay = imgTd.getTimeDelay(i)
        theta = imgTd.getImagePointPosition(imgNum, ptNum)
        print(theta/ANGLE_ARCSEC, timeDelay)
        beta = l.traceTheta(Ds, Dds, theta)
        betaAvg += beta
        thetas.append(theta)
        #print(theta/ANGLE_ARCSEC, beta/ANGLE_ARCSEC)

    betaAvg /= imgTd.getNumberOfTimeDelays()
    thetas = np.array(thetas)

    tds = l.getTimeDelay(zd, Ds, Dds, thetas, betaAvg)
    tds -= tds.min()
    # The function call returns SI units, wheras the input stored in the image info is in days
    print("")
    print("Time delays calculated from model")
    print(tds/(60*60*24))

print("Comparison for non-corrected lens:")
compareTimeDelays(baseLens)

print("Comparison for corrected lens:")
compareTimeDelays(correctedLens)

