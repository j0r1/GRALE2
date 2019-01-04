import grale.plotutil as plotutil
from grale.constants import *
import matplotlib.pyplot as plt

plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

import pickle
lensInfo = pickle.loads(open("lensinfo.bin", "rb").read())

plt.figure(figsize=(13,4))
plt.subplot(1,3,1)
plotutil.plotDensityContours(lensInfo, levels=20)
plt.gca().invert_xaxis()
plt.gca().set_title("Density contours")

plt.subplot(1,3,2)
plotutil.plotAverageDensityProfile(lensInfo, 50*ANGLE_ARCSEC)
plt.gca().set_title("Density profile")

plt.subplot(1,3,3)
plotutil.plotIntegratedMassProfile(lensInfo, 50*ANGLE_ARCSEC,
                                   massUnit = MASS_SUN)
plt.gca().set_title("Mass profile")

plt.show()

