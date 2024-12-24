from grale.all import *
import pickle

zd = 0.039
cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)

imgList = images.readInputImagesFile("images.txt", True)

w = 1*ANGLE_ARCSEC
initialParamDesc = {
    "type": "CompositeLens",
    "params": [
        {
            "x": { "initmin": -w, "initmax": +w, "hardmin": -w, "hardmax": +w },
            "y": { "initmin": -w, "initmax": +w, "hardmin": -w, "hardmax": +w },
            "factor": 1,
            "angle": { "initmin": 0, "initmax": 180, "hardmin": 0, "hardmax": 180 },
            "lens": {
                "type": "NSIELens",
                "params": {
                    "velocityDispersion": [ 200000 ],
                    "ellipticity": [ 0.8 ],
                    "coreRadius": 0.0001*ANGLE_ARCSEC, # Use small core radius to mimic a SIE
                }
            }
        }
    ]
}

def initialInversion(imgList, paramDesc):

    # This configures the inversion code to continuously add some random noise
    # to the input positions, to prevent overfitting
    imgListWithUncert = images.addPositionUncertainty(imgList, 0.05*ANGLE_ARCSEC)

    # Create the inversion workspace (the size doesn't really matter for a
    # parametric inversion, but is still a required parameter)
    iws = inversion.InversionWorkSpace(zd, 10*ANGLE_ARCSEC)
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
    fitParams = { "fitness_bayesweaklensing_stronglenssigma": 0.05*ANGLE_ARCSEC}

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
