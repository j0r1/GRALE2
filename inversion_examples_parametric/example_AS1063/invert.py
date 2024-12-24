from grale.all import *
import os
import pickle

# Lens redshift and cosmological model

zd = 0.3480
cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)
Dd = cosm.getAngularDiameterDistance(zd)

# Load the image information from "imawcs.mul"

def la(line):
    parts = line.split()
    src, img = map(int, parts[0].split("."))
    ra, dec = float(parts[1])*ANGLE_DEGREE, float(parts[2])*ANGLE_DEGREE
    z = float(parts[3])

    return { "srcnr": src, "imgnr": img, "x": ra, "y": dec, "z": z }

ctr = V(342.183210, -44.530878)*ANGLE_DEGREE
imgList = images.readInputImagesFile("./images.txt", True, la, ctr)

# This is the description of the lens model, some parameters are fixed,
# some are allowed to change and will be optimized.

initialParamDesc = {
    # A simple PIEMDLens always has a fixed location at (0,0), and fixed
    # orientation. To be able to change this, we'll use this inside a
    # CompositeLens model that can have multiple components (but we'll
    # use only one), re-centered and rotated.
    "type": "CompositeLens",
    "params": [
        {
            # A component's strength could be scaled by some factor, but
            # here we'll just set it to the fixed value of 1
            "factor": 1,
            # The center of the PIEMDLens below will be allowed to change
            # between these initial values (no hard bounds are set, so in
            # principle these values could go outside of these "initmin"
            # and "initmax" values)
            "x": { "initmin": -4*ANGLE_ARCSEC, "initmax": 4*ANGLE_ARCSEC },
            "y": { "initmin": -4*ANGLE_ARCSEC, "initmax": 4*ANGLE_ARCSEC },
            # Similar for the angle (in degrees), but here we do set hard
            # bounds
            "angle": { "initmin": 0, "initmax": 180, "hardmin": 0, "hardmax": 180 },
            # Here comes the actual lens model that's being recentered and
            # re-oriented. This is our PIEMDLens
            "lens": {
                "type": "PIEMDLens",
                "params": {
                    "centraldensity": { "initmin": 1, "initmax": 10 },
                    # Allow core radius to vary, but fix the scale radius
                    "coreradius": { "initmin": 1*ANGLE_ARCSEC, "initmax": 35*ANGLE_ARCSEC },
                    "scaleradius": 2000*ANGLE_ARCSEC,
                    "epsilon": { "initmin": 0.1, "initmax": 0.9 }
                }
            }
        }
    ]
}

# This is a helper function that will be called later. It performs the first
# inversion, not using MCMC but using the genetic algorithm to find a single
# best solution for this parametric description
def initialInversion(imgList, paramDesc):

    # This configures the inversion code to continuously add some random noise
    # to the input positions, to prevent overfitting
    imgListWithUncert = images.addPositionUncertainty(imgList, 0.05*ANGLE_ARCSEC)

    # Create the inversion workspace (the size doesn't really matter for a
    # parametric inversion, but is still a required parameter)
    iws = inversion.InversionWorkSpace(zd, 100*ANGLE_ARCSEC)
    for i in imgListWithUncert:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")

    # Invert!
    startLens, _, _, _, _ = iws.invertParametric(paramDesc, 512, eaType = "GA",
                                                 useImagePositionRandomization = True)
    startLens.save("startInv.lensdata")

    # This routine will create a modified parametric description, where
    # the initial parameter ranges will be close to the value from the
    # solution that was found. This can then be used in an MCMC run to
    # explore parameter uncertainties
    return paramdesc.refineParametricDescription(paramDesc, startLens, 0.001)

# This is the second helper function, that will be called with the refined
# parametric description. It uses the Goodman-Weare algorithm (same as 'emcee')
# to explore the parameter space
def mcmcInversion(imgList, paramDesc):

    iws = inversion.InversionWorkSpace(zd, 100*ANGLE_ARCSEC)
    for i in imgList:
        # Here we want to use the logprob based on the deviations in the
        # image plane. This is done by specifying that the images are used
        # in the bayesian code, for the strong lensing contribution.
        iws.addImageDataToList(i["imgdata"], i["z"], "bayesstronglensing")

    # Well run the algorithm for 2000 generation, half of which will be
    # ignored as burn-in samples
    totalGenerations = 5000
    burnInGenerations = totalGenerations//2

    # We'll also specify the file to which the samples are written
    mcmcParams = { "burningenerations": burnInGenerations, "samplesfilename": "samples.dat" }

    # This specifies that all strong lensing images are considered to have an
    # uncertainty of 1 arcsec
    fitParams = { "fitness_bayesweaklensing_stronglenssigma": 1*ANGLE_ARCSEC}

    # Run the MCMC algorithm! It will move 512 'walkers' around, for 2000 generations.
    lens, _, _, reducedParamInfo, _ = iws.invertParametric(paramDesc, 512, maximumGenerations = totalGenerations, eaType = "MCMC", 
                                                           geneticAlgorithmParameters=mcmcParams, fitnessObjectParameters=fitParams)

    lens.save("mcmcInv.lensdata")

    # To interpret the samples that are written to file, we also write out the
    # information about the parameters that are being varied. These parameters
    # typically use some find of scale factor, to allow for 32-bit floating point
    # calculations, and it's those scaled values that are written to the file.
    # The 'paramInfo' list contains the scale factors to convert the values from
    # file to the real values.
    pickle.dump(reducedParamInfo, open("paraminfo.dat", "wb"))

# Run the initial inversion to get a good starting point for the MCMC part
newParamDesc = initialInversion(imgList, initialParamDesc)
# Do the actual MCMC inversion
mcmcInversion(imgList, newParamDesc)


