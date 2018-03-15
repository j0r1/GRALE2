import grale.plotutil as plotutil
from grale.constants import *
import matplotlib.pyplot as plt

plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

import pickle
lensInfo = pickle.loads(open("lensinfo.bin", "rb").read())

import grale.images as images
import numpy as np

V = lambda x,y: np.array([x,y], dtype=np.double)

src = images.CircularSource(V(5.6,-5)*ANGLE_ARCSEC, 1.1*ANGLE_ARCSEC, fade=True)

plt.figure(figsize=(5,5))
plotutil.plotImagePlane(lensInfo, [src])
plt.gca().invert_xaxis()
plt.show()
