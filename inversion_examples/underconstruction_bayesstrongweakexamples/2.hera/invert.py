# This is very similar to the previous SL/WL example, but now using
# the Hera data. Many things were not known, and can certainly be
# improved.

from grale.all import *
import os

z_lens = 0.5073
cosm = cosmology.Cosmology(0.72, 0.24, 0, 0.76)
cosmology.setDefaultCosmology(cosm)
D = cosm.getAngularDiameterDistance

# Helper function to read the WL data from the catalog
def getShearInfo():
    
    cx = 4500
    cy = 4500
    arcsecPerPixel = 0.2
    
    imgDat = images.ImagesData(1, shear=True, shearsigma=True, redshift=True, redshiftsigma=True)
    
    for l in open("finalshear_z.cat"):
        if l.strip().startswith("#"):
            continue
        xpix, ypix, gx, gy, z = map(float, l.split())
        
        x = (xpix-cx)*arcsecPerPixel*ANGLE_ARCSEC
        y = (ypix-cy)*arcsecPerPixel*ANGLE_ARCSEC
        
        # Note: based on how the reconstruction looks, it seems that using -x, -y
        # should be used. This makes sense for the x direction if the RA axis is
        # used (which points left). And for the y direction if pixel y coordinates
        # go from top to bottom (as is quite common in pixel-based content)
        imgDat.addPoint(0, [-x, -y], shear=[gx, gy], shearsigma=[0,0], # TODO: add uncertainty
                        redshift=z, redshiftsigma=0.5) # TODO: what's a good redshift sigma?
    
    return imgDat

useTestRun = False
popSize = 16 if useTestRun else 512

targetMassScale = "auto"

useRMS = True

scaleSearchType = "wide"
ignoreWLMass = True if scaleSearchType == "regular" else False

weakSubDiv = 24
weakMassScale = 1e15*MASS_SUN

sheetType = "nosheet"

inversion.setDefaultInverter("threads")

weakRegSize = 30*ANGLE_ARCMIN 
strongRegSize = 120*ANGLE_ARCSEC
iws = inversion.InversionWorkSpace(z_lens, strongRegSize)

def lineAnalyzer(line):
    x, y, src, img, z = map(float, line.split())
    return { "x": x*ANGLE_ARCSEC, "y": y*ANGLE_ARCSEC, "z": z, "srcnr": src, "imgnr": img }

imgList = images.readInputImagesFile("multimages.txt", True, lineAnalyzer) 
for i in imgList:
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
    if useRMS:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointgroupimages")

if weakSubDiv > 0:
    iws.addImageDataToList(getShearInfo(), None, "bayesellipticities")

def baDist(N, low, high):
    xs = np.linspace(0, 1, N)
    ys = np.array([1.0 if x > low and x < high else (x/low if x < low else (1-x)/(1.0-high)) for x in xs ])
    return ys

fop = {
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

