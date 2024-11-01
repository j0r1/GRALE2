from grale.all import *
import matplotlib.pyplot as plt
import numpy as np
import pickle

feedback.setDefaultFeedback("none")

np.random.seed(12345)

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)
D = cosm.getAngularDiameterDistance

zd = 0.4
Dd = D(zd)

trueLens = lenses.CompositeLens(Dd, [
    {
        "lens": lenses.NSIELens(Dd, { "velocityDispersion": 700000, "ellipticity": 0.8, "coreRadius": 0.5*ANGLE_ARCSEC }),
        "x": -10*ANGLE_ARCSEC, "y": 0, "angle": 30, "factor": 1,
    },
    {
        "lens": lenses.NSIELens(Dd, { "velocityDispersion": 900000, "ellipticity": 0.75, "coreRadius": 1.5*ANGLE_ARCSEC }),
        "x": 5*ANGLE_ARCSEC, "y": 10*ANGLE_ARCSEC, "angle": -60, "factor": 1,
    }
])

trueLens.save("truelens.lensdata");

li = plotutil.LensInfo(trueLens, size=80*ANGLE_ARCSEC, zd=zd, zs=2)

sigma = 0.1*ANGLE_ARCSEC

numSources = 20
srcPos = []
imgListNoNoise = []
imgListNoise = []
while len(srcPos) < numSources:
    x = np.random.uniform(-15, 15)*ANGLE_ARCSEC
    y = np.random.uniform(-15, 15)*ANGLE_ARCSEC
    z = np.random.uniform(zd*1.5, 4.5)
    li.setSourceRedshift(z)
    ip = li.getImagePlane()
    imgs = ip.traceBeta(V(x,y))
    if len(imgs) > 1:
        imgDatNoise = images.ImagesData(len(imgs))
        imgDatNoNoise = images.ImagesData(len(imgs))
        for imgNum, pt in enumerate(imgs):
            imgDatNoNoise.addPoint(imgNum, pt)

            pt = np.array(pt) + np.random.normal(loc=0, scale=sigma, size=2) # Add some small noise
            imgDatNoise.addPoint(imgNum, pt)

        imgListNoise.append({"imgdata": imgDatNoise, "z": z})
        imgListNoNoise.append({"imgdata": imgDatNoNoise, "z": z})
        srcPos.append(V(x,y))
        print(f"{len(srcPos)}/{numSources}")

pickle.dump(imgListNoise, open("imglist_noise.pickle", "wb"))
pickle.dump(imgListNoNoise, open("imglist_nonoise.pickle", "wb"))
pickle.dump(srcPos, open("srcpos.pickle", "wb"))

plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
plotutil.plotImagePlane(li)
plt.subplot(2,2,2)
plotutil.plotDensityContours(li, levels=np.arange(0.5, 15, 0.5))
plt.subplot(2,2,3)
plotutil.plotImagesData(imgListNoNoise)
plt.show()
