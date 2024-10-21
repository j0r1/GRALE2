from grale.all import *
import pickle
import pprint

inverters.debugOutput = True
inverters.debugDirectStderr = True

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)
D = cosm.getAngularDiameterDistance

zd = 0.4
Dd = D(zd)

imgList = pickle.load(open("imglist.pickle", "rb"))
for i in imgList:
    i["images"] = i["imgdata"] # Internally another key is used
    i["Ds"] = D(i["z"])
    i["Dds"] = D(zd, i["z"])
    i["params"] =  {"type":"pointimages"}

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
result = inversion.invertParametric(imgList, lensDescription, zd, Dd, 128, inverter="mpi:2")
pprint.pprint(result)
