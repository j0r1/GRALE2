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

# The entry in this array will be used to accumulate the pixels of the 
# different image planes in. Here, the numpy array is wrapped in an array
# to avoid a 'local variable is referenced before assignment' error
sumPlane = [ None ]

# This function will intercept the pixels before they're plotted using
# matplotlib's 'imshow'. What is returned here will be what is actually
# shown, which is the sum of the individual image planes.
def f(x):
    if sumPlane[0] is None: sumPlane[0] = np.zeros(x.shape)
    sumPlane[0] += x
    return sumPlane[0]

# Here, we tell the plotImagePlane function to call 'f' right before
# plotting the image, thereby allowing the pixels to be accumulated
plt.figure(figsize=(5,5))
plotutil.plotImagePlane(lensInfo, [src], processRenderPixels=f)

# For the second source, we still need to set different source
# distances. We'll also plot the caustics and critical lines in a
# slightly different color.
lensInfo.setSourceDistances(Ds2, Dds2)
plotutil.plotImagePlane(lensInfo, [src2], processRenderPixels=f,
                        caustColor="darkblue", critColor="darkred");

plt.gca().invert_xaxis()
plt.show()

