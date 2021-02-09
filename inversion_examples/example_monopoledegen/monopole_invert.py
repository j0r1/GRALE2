from grale.all import *
import random

inversion.setDefaultInverter("mpi")

# These are the lens parameters that were used in making the model
cosm = cosmology.Cosmology(0.71, 0.27, 0, 0.73)
cosmology.setDefaultCosmology(cosm)

z_lens = 0.395
Dd = cosm.getAngularDiameterDistance(z_lens)

# To make experimenting with this script easier, an approximate lens
# model is used. The full model can be used of course, and while the
# inversion/optimization itself will have the same speed, initializing
# it will take longer. The approximation is based on the deflection
# angles on a 511x511 grid covering a 2 arcmin region.
baseLens = lenses.GravitationalLens.load("approxlens.lensdata")

srcs = [ images.ImagesData.load("source1entire.imgdata"),
         images.ImagesData.load("source2entire.imgdata") ]

# The createMonopoleBasisFunctions function below creates the monopole 
# basis functions of which the weights will need to be optimized. Using
# a callback for each cell, we'll also create the ImagesData instance that
# will be used for the numeric estimation of the gradient.

gradImg = images.ImagesData(3)
    
def cellCenterCallback(ctr, cellSize, hasBasisFunction, gradImg):
    cellPartForNumericGradient = 1000.0
    dx = cellSize/cellPartForNumericGradient

    gradImg.addPoint(0, ctr)
    gradImg.addPoint(1, ctr + V(dx, 0))
    gradImg.addPoint(2, ctr + V(0, dx))
        
# These are the same settings as used for the figure in the CL0024 article
basisFunctions = util.createMonopoleBasisFunctions(srcs, Dd, 31, 1.2*ANGLE_ARCMIN,
                                              cellCenterCallback=cellCenterCallback,
                                              cellCenterCallbackState=gradImg)

# This is not really relevant for the kind of optimization routine we'll be
# using
totalFakeMass = sum([ b["mass"] for b in basisFunctions ])

iws = inversion.InversionWorkSpace(z_lens, 1.6*ANGLE_ARCMIN)
iws.addImageDataToList(gradImg, 1.0, "kappagradientgrids")
iws.addBasisFunctions(basisFunctions)

# Originally, in 2008, a separate genetic algorithm was used, but thanks to
# some rewrites in the inversion code, the same algorithm can be used for this
# type of optimization as well.
# The original code always used 2000 generations, which is what we'll do here
# as well. The maximumGenerations is of course relevant, but 
# 'convergence_history_size' as well. To check if the GA should be stopped,
# the current fitness is compared to that a number of generations ago. By
# setting this number to 2000 as well, we make sure that the algorithm can't
# stop earlier than after 2000 generations
popSize = 512
numGen = 2000

lens, _, _ = iws.invertBasisFunctions(popSize,
    massScaleSearchType="nosearch",
    massScale=totalFakeMass, # Not used for 'nosearch' but something must be set
    fitnessObjectParameters={
        "convergence_history_size": numGen,
    },
    maximumGenerations=numGen,
    baseLens=baseLens,
    allowNegativeValues=True
    )

lens.save("monopoles.lensdata")

