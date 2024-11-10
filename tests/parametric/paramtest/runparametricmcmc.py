from grale.all import *
import pickle
import pprint

inverters.debugOutput = True
inverters.debugDirectStderr = True

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)

zd = 0.4

lensDescription = {
    "type": "CompositeLens",
    "params": [
        { 
            "x": [-10*ANGLE_ARCSEC,0.2,0.2 ], "y": [1*ANGLE_ARCSEC, 0.5, 3], "factor": 1, "angle": [30, 0.5],
            "lens": {
                #"type": "NSIELens",
                #"params": {
                #    "velocityDispersion": [700000, 0.2, 0.5],
                #    "ellipticity": [0.8, 0.1],
                #    "coreRadius": [0.5*ANGLE_ARCSEC, 0.5, 0.9]
                #}
                "type": "PIEMDLens",
                "params": {
                    "epsilon": { "initmin": 0.1, "initmax": 0.9 },
                    "coreradius": [ 1*ANGLE_ARCSEC ],
                    "scaleradius": [ 20*ANGLE_ARCSEC ],
                    "centraldensity": { "initmin": 1, "initmax": 50 },
                }
            }
        },
        {
            "x": [5*ANGLE_ARCSEC,0.2,0.2], "y": [10*ANGLE_ARCSEC, 0.2, 0.2], "factor": 1, "angle": [-60, 0.5],
            "lens": {
                #"type": "NSIELens",
                #"params": {
                #    "velocityDispersion": [900000, 0.2, 0.5],
                #    "ellipticity": [0.75, 0.1],
                #    "coreRadius": [1.5*ANGLE_ARCSEC, 0.5, 0.9]
                #}
                "type": "PIEMDLens",
                "params": {
                    "epsilon": { "initmin": 0.1, "initmax": 0.9 },
                    "coreradius": [ 1*ANGLE_ARCSEC ],
                    "scaleradius": [ 20*ANGLE_ARCSEC ],
                    "centraldensity": { "initmin": 1, "initmax": 50 },
                }
            }
        }
    ]
}

#inverters.debugCaptureProcessCommandsFile = "dumpcommunication.dat"
#inversion.setDefaultInverter("threads:8")
#import os
#os.environ["GRALE_DEBUG_SEED"] = "12345"

for imgListFn, solFn, sampFn, paramFn in [ ("imglist_noise.pickle", "sol_mcmc_noise.lensdata", "samples_noise.dat", "params_noise.dat"),
                                           ("imglist_nonoise.pickle", "sol_mcmc_nonoise.lensdata", "samples_nonoise.dat", "params_nonoise.dat") ]:

    imgList = pickle.load(open(imgListFn, "rb"))

    # Do a non-MCMC inversion to get a good starting point as lens
    iws = inversion.InversionWorkSpace(zd, 10*ANGLE_ARCSEC)
    for i in imgList:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")

    lens, _, _, paramInfo = iws.invertParametric(lensDescription, 512)
    usedParameters = [ x["name"] for x in paramInfo ]

    def cvf(value, lensNameList, paramName, uniqueParamName, allParams):
        if uniqueParamName in usedParameters:
            print("Using", uniqueParamName)
            return [ value, 0.01 ] # Start in the neighbourhood of the current best parameters
        return value

    refinedLensDescription = paramdesc.createParametricDescription(lens, convertValueFunction=cvf)
    refinedLensDescription = eval(refinedLensDescription)

    iws = inversion.InversionWorkSpace(zd, 10*ANGLE_ARCSEC)
    for i in imgList:
        iws.addImageDataToList(i["imgdata"], i["z"], "bayesstronglensing")

    lens, _, _, paramInfo = iws.invertParametric(refinedLensDescription, 512, maximumGenerations=2000, eaType = "MCMC",
                              fitnessObjectParameters = { "fitness_bayesweaklensing_stronglenssigma": 0.05*ANGLE_ARCSEC},
                              geneticAlgorithmParameters = { "annealgenerationsscale": 900, "burningenerations": 1000,
                                                             "samplesfilename": sampFn})
    lens.save(solFn)
    pickle.dump(paramInfo, open(paramFn, "wb"))

