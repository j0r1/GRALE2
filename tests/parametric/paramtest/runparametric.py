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
               # "type": "NSIELens",
               # "params": {
               #     "velocityDispersion": [700000, 0.2, 0.5],
               #     "ellipticity": [0.8, 0.1],
               #     "coreRadius": [0.5*ANGLE_ARCSEC, 0.5, 0.9]
               # }
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
#null = images.createGridTriangles(V(-100,-100)*ANGLE_ARCSEC, V(100,100)*ANGLE_ARCSEC, 48,48)

#import os
#os.environ["GRALE_DEBUG_SEED"] = "12345"

for imgListFn, solFn in [ ("imglist_noise.pickle", "sol_jade_noise.lensdata"),
                          ("imglist_nonoise.pickle", "sol_jade_nonoise.lensdata") ]:

    imgList = pickle.load(open(imgListFn, "rb"))

    iws = inversion.InversionWorkSpace(zd, 10*ANGLE_ARCSEC)
    for i in imgList:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
        #iws.addImageDataToList(null, i["z"], "pointnullgrid")
        #iws.addImageDataToList(i["imgdata"], i["z"], "pointgroupimages")

    lens, _, _, _, _ = iws.invertParametric(lensDescription, 512)
    lens.save(solFn)
