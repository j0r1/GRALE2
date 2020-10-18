from grale.all import *

plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

z_lens = 0.4
iws = inversion.InversionWorkSpace(z_lens, 250*ANGLE_ARCSEC, cosmology=cosmology.Cosmology(0.7, 0.27, 0, 0.73))

imgList = images.readInputImagesFile("ex1_inputpoints.txt", True)
for i in imgList:
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")

lens1 = lenses.GravitationalLens.load("ex1_inv1.lensdata")
lens2 = lenses.GravitationalLens.load("ex1_inv2.lensdata")

bpImgs1 = iws.backProject(lens1)
bpImgs2 = iws.backProject(lens2)

import matplotlib.pyplot as plt

plt.figure(figsize=(15,5))
plt.subplot(1,3,1)
plt.title("Input point images")
plotutil.plotImagesData(imgList)
plt.gca().set_aspect("equal")

plt.subplot(1,3,2)
plt.title("Back-projected with 'lens1' (uniform grid)")
plotutil.plotImagesData(bpImgs1)
plt.gca().set_aspect("equal")

plt.subplot(1,3,3)
plt.title("With 'lens2' (subdivision grid)")
plotutil.plotImagesData(bpImgs2)
plt.gca().set_aspect("equal")

plt.tight_layout()
plt.show()
