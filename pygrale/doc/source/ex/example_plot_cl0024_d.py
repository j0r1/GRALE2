import grale.plotutil as plotutil
import grale.cosmology as cosmology
from grale.constants import *
import matplotlib.pyplot as plt

cosm = cosmology.Cosmology(0.71, 0.27, 0, 0.73)
zd, zs = 0.395, 1.675
Ds = cosm.getAngularDiameterDistance(zs)
Dds = cosm.getAngularDiameterDistance(zd, zs)

plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

import pickle
lensInfo = pickle.loads(open("lensinfo.bin", "rb").read())

import grale.images as images
import numpy as np

V = lambda x,y: np.array([x,y], dtype=np.double)
src = images.CircularSource(V(5.6,-5)*ANGLE_ARCSEC, 1.1*ANGLE_ARCSEC, fade=True)

# Get distances and shape for a second source
zs2 = 1.3
Ds2 = cosm.getAngularDiameterDistance(zs2)
Dds2 = cosm.getAngularDiameterDistance(zd, zs2)
src2 = images.CircularSource(V(0.21,-9.26)*ANGLE_ARCSEC, 0.75*ANGLE_ARCSEC, fade=True)

# Here, we tell the plotImagePlane function to call 'f' right before
# plotting the image, thereby allowing the pixels to be accumulated
plt.figure(figsize=(5,5))
plotutil.plotImagePlane(lensInfo, [{"shape": src, "Ds": Ds, "Dds": Dds},
                                   {"shape": src2, "Ds": Ds2, "Dds": Dds2}],
                        critColor = [ 'red', 'darkred'],
                        caustColor = [ 'blue', 'darkblue'])

plt.gca().invert_xaxis()
plt.show()

