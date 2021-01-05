from grale.all import *
import os
import glob
import pprint
# This helper file contains a function that provides the WL data
from weakdata import getWeakData

import random
print("RNG State:")
print(random.getstate())

z_lens = 0.4
cosm = cosmology.Cosmology(0.7, 0.27, 0, 0.73)
cosmology.setDefaultCosmology(cosm)

strongSize = 200*ANGLE_ARCSEC
weakRegSize = 30*ANGLE_ARCMIN 

# This is for testing, to see what happens if the weak region size is larger
doubleSize = False
if doubleSize:
    weakRegSize *= 2

weakMassScale = 1e15*MASS_SUN

useTestRun = False
popSize = 16 if useTestRun else 512

targetMassScale = "auto"

# Adds an extra fitness measure, I'm not entirely sure how much difference 
# this makes, but it does seem like a small improvement
useRMS = True

# Is passed to the getWeakData function, if set it removes parts of the
# data (run the weakdata.py script to see which parts)
partialWL = False

scaleSearchType = "wide" # or "regular"
# If a (narrow) mass scale search based on only the SL data is used,
# we won't count the WL mass.

# While at first I thought that ignoring mass in the WL region would
# be a good idea, this has the side effect that there are no bounds
# on how much mass there can be in the entire system. I now generally
# tend to use this "wide" approach.
ignoreWLMass = True if scaleSearchType == "regular" else False

weakSubDiv = 24 # Set to 0 to disable WL
ellImg = getWeakData(partial = partialWL) if weakSubDiv else None

sheetType = "nosheet" # or "genome"

inversion.setDefaultInverter("mpi")

iws = inversion.InversionWorkSpace(z_lens, strongSize)
for i in images.readInputImagesFile("images.txt", True):
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
    if useRMS:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointgroupimages")

if ellImg:
    iws.addImageDataToList(ellImg, None, "bayesellipticities")

# This provides roughly the same uniform distribution as used when simulating
# the WL data. Instead of being 0 below 0.1 and above 0.9 however, this goes
# to zero linearly. This seems to make the inversion much better behaved.
def baDist(N, low, high):
    xs = np.linspace(0, 1, N)
    ys = np.array([1.0 if x > low and x < high else (x/low if x < low else (1-x)/(1.0-high)) for x in xs ])
    return ys

# This distribution will be passed to the bayes WL fitness using fitnessObjectParameters
# below.
fop = { 
    "fitness_bayesweaklensing_b_over_a_distribution": baDist(1024, 0.1, 0.9),
}

# For a quick test to see if all code executes
if useTestRun:
    print("INFO: USING TEST RUN")
    iws.setDefaultInversionArguments(maximumGenerations = 2)

prevLens, subDivStart = None, 100
bestLens, bestFitness = None, None
# These are the same 5 steps as in all inversion examples, but now
# using a for loop
for idx in range(1,6):

    strongSubDivInf = 15 if prevLens is None else (prevLens, subDivStart, subDivStart+100)

    # This is a convenience function, it combines strong and weak grid
    # basis functions. The docs provide more information.
    iws.setStrongAndWeakBasisFunctions(
        strongSubDivInf, weakSubDiv, weakRegSize, 
        weakMassScale=weakMassScale, ignoreWLMassInMassScaleSearch=ignoreWLMass)

    lens, fitness, fitdesc = iws.invertBasisFunctions(popSize,
                                                      sheetSearch = sheetType, fitnessObjectParameters = fop,
                                                      massScaleSearchType = scaleSearchType,
                                                      massScale = targetMassScale)
    
    fileName = f"inv{idx}.lensdata" 
    lens.save(fileName)
    print(f"LENSFITNESS: {fileName} {fitness}")

    prevLens = lens
    subDivStart += 200

    if bestLens is None or fitness < bestFitness:
        bestLens = lens
        bestFitness = fitness

print(f"BESTFITNESS: {bestFitness}")
bestLens.save("best.lensdata")

