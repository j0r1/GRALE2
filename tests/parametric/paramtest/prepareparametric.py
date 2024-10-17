from grale.all import *
import pprint

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)
D = cosm.getAngularDiameterDistance

zd = 0.4
Dd = D(zd)

lensDescription = {
    "type": "CompositeLens",
    "params": [
        { 
            "x": [-9.5*ANGLE_ARCSEC ], "y": [1*ANGLE_ARCSEC], "factor": 1, "angle": [20],
            "lens": {
                "type": "NSIELens",
                "params": {
                    "velocityDispersion": [ 700000 ],
                    "ellipticity": [0.85],
                    "coreRadius": [1*ANGLE_ARCSEC]
                }
            }
        },
        {
            "x": [5.5*ANGLE_ARCSEC ], "y": [11*ANGLE_ARCSEC], "factor": 1, "angle": [-50],
            "lens": {
                "type": "NSIELens",
                "params": {
                    "velocityDispersion": [ 1000000 ],
                    "ellipticity": [0.9],
                    "coreRadius": [1*ANGLE_ARCSEC]
                }
            }
        }
    ]
}

desc = paramdesc.analyzeParametricLensDescription(lensDescription, Dd, 0.1)
pprint.pprint(desc)

templateParams = desc["floatparams"]
varParams = desc["variablefloatparams"]
offsets = [ x["offset"] for x in varParams ]
initMin = [ x["initialrange"][0] for x in varParams ]
initMax = [ x["initialrange"][1] for x in varParams ]
hardMin = [ x["hardlimits"][0] for x in varParams ]
hardMax = [ x["hardlimits"][1] for x in varParams ]

desc["templatelens"].save("templatelens.lensdata")
with open("inversionparams.txt", "wt") as f:
    f.write(f"{desc["scales"]["deflectionscale"]:.15g} {desc["scales"]["potentialscale"]:.15g}\n")
    f.write(f"{len(offsets)}\n")
    f.write(" ".join(list(map(str, offsets))) + "\n")
    f.write(" ".join(list(map(lambda x:f"{x:.8g}", initMin))) + "\n")
    f.write(" ".join(list(map(lambda x:f"{x:.8g}", initMax))) + "\n")
    f.write(" ".join(list(map(lambda x:f"{x:.8g}", hardMin))) + "\n")
    f.write(" ".join(list(map(lambda x:f"{x:.8g}", hardMax))) + "\n")
