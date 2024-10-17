from grale.all import *
import matplotlib.pyplot as plt
import numpy as np

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

li = plotutil.LensInfo(trueLens, size=80*ANGLE_ARCSEC, zd=zd, zs=2)

if 1:
    numSources = 20
    srcPos = []
    imgList = []
    while len(srcPos) < numSources:
        x = np.random.uniform(-15, 15)*ANGLE_ARCSEC
        y = np.random.uniform(-15, 15)*ANGLE_ARCSEC
        z = np.random.uniform(zd*1.5, 4.5)
        li.setSourceRedshift(z)
        ip = li.getImagePlane()
        imgs = ip.traceBeta(V(x,y))
        if len(imgs) > 1:
            imgDat = images.ImagesData(len(imgs))
            for imgNum, pt in enumerate(imgs):
                imgDat.addPoint(imgNum, pt)

            imgList.append({"imgdata": imgDat, "z": z})
            srcPos.append(V(x,y))
            print(f"{len(srcPos)}/{numSources}")

with open("images.txt", "wt") as f:
    for i in imgList:
        z = i["z"]
        imgDat = i["imgdata"]
        for imgNum in range(imgDat.getNumberOfImages()):
            pt = imgDat.getImagePointPosition(imgNum, 0)/ANGLE_ARCSEC
            f.write(f"{pt[0]:g} {pt[1]:g} {z:g}\n")
        f.write("\n")


plt.figure(figsize=(10,10))
plt.subplot(2,2,1)
plotutil.plotImagePlane(li)
plt.subplot(2,2,2)
plotutil.plotDensityContours(li, levels=np.arange(0.5, 15, 0.5))
plt.subplot(2,2,3)
plotutil.plotImagesData(imgList)
plt.show()