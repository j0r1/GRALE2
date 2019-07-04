from example_lens_dl import *
import grale.plotutil as plotutil
import grale.renderers as renderers
import grale.feedback as feedback
import grale.cosmology as cosmology
from grale.constants import *
import matplotlib.pyplot as plt

# Let's create an alias for this class
LI = plotutil.LensInfo

# Set defaults
renderers.setDefaultMassRenderer("threads")
renderers.setDefaultLensPlaneRenderer("threads")
feedback.setDefaultFeedback("stdout")
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

# For the generation of this model of CL0024, the
# following cosmology and redshifts were used
cosm = cosmology.Cosmology(0.71, 0.27, 0, 0.73)
zd, zs = 0.395, 1.675
Ds = cosm.getAngularDiameterDistance(zs)
Dds = cosm.getAngularDiameterDistance(zd, zs)

# Inform the LensInfo object that we're interested in a region that's
# 120" wide, with the specified angular diameter distances.
lensInfo = LI(cl0024Model, size=120*ANGLE_ARCSEC, Ds=Ds, Dds=Dds)

try:
    # This is actually mainly for convenience in generating this
    # documentation: at the end of the file we'll save the LensInfo
    # object, which by then contains calculated lens mapping and mass
    # density data, to a file. If this file already exists, we can
    # just use that to obtain the previously calculated data.
    #
    # This makes generating the documentation much less cumbersome,
    # as calculating the lens mapping and mass density map takes several
    # minutes.
    import pickle
    lensInfo = pickle.loads(open("lensinfo.bin", "rb").read())
except Exception as e:
    print(e)


plt.figure(figsize=(12,4))

# Create a plot of the image plane for the cl0024 model. Since the
# right ascension axis goes from right to left, we'll invert the
# X-axis direction.
plt.subplot(1,3,1)
plotutil.plotImagePlane(lensInfo)
plt.gca().invert_xaxis()

# And a plot for the density
plt.subplot(1,3,2)
plotutil.plotDensity(lensInfo)
plt.gca().invert_xaxis()

# Let's also combine the two plots
plt.subplot(1,3,3)
plotutil.plotDensity(lensInfo, cmap="gray_r")
plotutil.plotImagePlane(lensInfo, bgRgb=(0,0,0,0))
plt.gca().invert_xaxis()

plt.show()

# As mentioned before, here we save the generated LensInfo object, to a
# file.
open("lensinfo.bin", "wb").write(pickle.dumps(lensInfo))

