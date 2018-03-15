from example_lens import *
import grale.plotutil as plotutil
import matplotlib.pyplot as plt

plt.figure(figsize=(10,5))
# Set the unit in which to _display_ plots to one arcsec
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

# Create a plot of the image plane for this combined lens
plt.subplot(1,2,1)
plotutil.plotImagePlane(combLens)

# And a plot for the density
plt.subplot(1,2,2)

# This callback will be used to add a color bar
def addBar(a):
    plt.colorbar(a, fraction=0.046, pad=0.04)

plotutil.plotDensity(combLens, axImgCallback=addBar, vmax=10)
plt.show()
