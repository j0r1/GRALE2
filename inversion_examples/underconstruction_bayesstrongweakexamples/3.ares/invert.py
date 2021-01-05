# Similar to the previous two examples, using the Ares data
from grale.all import *
import os

import random
print("RNG State:")
print(random.getstate())

cosm = cosmology.Cosmology(0.704, 0.272, 0, 0.728)
cosmology.setDefaultCosmology(cosm)

z_lens = 0.5
weakRegSize = 30*ANGLE_ARCMIN 
strongRegSize = 170*ANGLE_ARCSEC
iws = inversion.InversionWorkSpace(z_lens, strongRegSize)

useTestRun = False
popSize = 16 if useTestRun else 512

targetMassScale = "auto"

useRMS = True

scaleSearchType = "wide"
ignoreWLMass = True if scaleSearchType == "regular" else False

weakSubDiv = 24
weakMassScale = 1e15*MASS_SUN

def getEllInfo():
    imgDat = images.ImagesData(1, shear=True, shearsigma=True, redshift=True, redshiftsigma=True)
    cx = 4500
    cy = 4500
    arcsecPerPixel = 0.2
    for l in open("finalshear.cat"):
        xpix, ypix, gx, gy = map(float, l.split())
        x = (xpix-cx)*arcsecPerPixel*ANGLE_ARCSEC
        y = (ypix-cy)*arcsecPerPixel*ANGLE_ARCSEC

        # Here we do the same -x, -y operation as for the Hera data. It doesn't
        # cause a clear difference though, the mass density is quite symmetric.
        imgDat.addPoint(0, [-x, -y], 
                        shear=[gx, gy], shearsigma=[0,0], # TODO: add uncertainties for ellipticities
                        redshift=0, redshiftsigma=0, # Completely unknown redshift
                        )

    return imgDat

sheetType = "genome"
inversion.setDefaultInverter("mpi")

def lineAnalyzer(line):
    x, y, src, img, z = map(float, line.split())
    return { "x": x*ANGLE_ARCSEC, "y": y*ANGLE_ARCSEC, "z": z, "srcnr": src, "imgnr": img }

imgList = images.readInputImagesFile("multimages.txt", True, lineAnalyzer) 
for i in imgList:
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
    if useRMS:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointgroupimages")

if weakSubDiv > 0:
    iws.addImageDataToList(getEllInfo(), None, "bayesellipticities")

def baDist(N, low, high):
    xs = np.linspace(0, 1, N)
    ys = np.array([1.0 if x > low and x < high else (x/low if x < low else (1-x)/(1.0-high)) for x in xs ])
    return ys

fop = {
    # TODO: need redshift distribution, what's good here? This is a simple uniform dist from 0 to 5
    "fitness_bayesweaklensing_zdist_values": np.array([1.0,1.0]),
    "fitness_bayesweaklensing_zdist_range": np.array([0.0,5.0]),
    # TODO: what's good here?
    "fitness_bayesweaklensing_b_over_a_distribution": baDist(1024, 0.1, 0.9),
}

# For a quick test to see if all code executes
if useTestRun:
    print("INFO: USING TEST RUN")
    iws.setDefaultInversionArguments(maximumGenerations = 2)

prevLens, subDivStart = None, 100
bestLens, bestFitness = None, None

for idx in range(1,6):

    strongSubDivInf = 15 if prevLens is None else (prevLens, subDivStart, subDivStart+100)
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
