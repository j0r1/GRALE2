from grale.all import *
import pickle
import pprint

inverters.debugOutput = True
inverters.debugDirectStderr = True

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)

zd = 0.4
imgList = pickle.load(open("imglist.pickle", "rb"))

lensDescription = {
    "type": "CompositeLens",
    "params": [
        { 
            "x": [-10*ANGLE_ARCSEC,0.2,0.2 ], "y": [1*ANGLE_ARCSEC, 0.5, 3], "factor": 1, "angle": [30, 0.5],
            "lens": {
                "type": "NSIELens",
                "params": {
                    "velocityDispersion": [700000, 0.2, 0.5],
                    "ellipticity": [0.8, 0.1],
                    "coreRadius": [0.5*ANGLE_ARCSEC, 0.5, 0.9]
                }
            }
        },
        {
            "x": [5*ANGLE_ARCSEC,0.2,0.2], "y": [10*ANGLE_ARCSEC, 0.2, 0.2], "factor": 1, "angle": [-60, 0.5],
            "lens": {
                "type": "NSIELens",
                "params": {
                    "velocityDispersion": [900000, 0.2, 0.5],
                    "ellipticity": [0.75, 0.1],
                    "coreRadius": [1.5*ANGLE_ARCSEC, 0.5, 0.9]
                }
            }
        }
    ]
}

#inverters.debugCaptureProcessCommandsFile = "dumpcommunication.dat"

inversion.setDefaultInverter("threads:8")

null = images.createGridTriangles(V(-100,-100)*ANGLE_ARCSEC, V(100,100)*ANGLE_ARCSEC, 48,48)

iws = inversion.InversionWorkSpace(zd, 10*ANGLE_ARCSEC)
for i in imgList:
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
    iws.addImageDataToList(null, i["z"], "pointnullgrid")
    #iws.addImageDataToList(i["imgdata"], i["z"], "pointgroupimages")

import os
os.environ["GRALE_DEBUG_SEED"] = "12345"

result = iws.invertParametric(lensDescription, 64, maximumGenerations=1000)
pprint.pprint(result)
