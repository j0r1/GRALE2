from example_lens_dl import *
import grale.plotutil as plotutil
import grale.renderers as renderers
from grale.constants import *
import numpy as np
import matplotlib.pyplot as plt
import pickle

# Use the available cores for calculating the mass densities, and set the
# plot units to arc seconds
renderers.setDefaultMassRenderer("threads")
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

# Create some aliases for these class names
LI = plotutil.LensInfo
DI = plotutil.DensInfo

# The model itself is a CompositeLens, and the parameters will contain
# the individual lens models
params = cl0024Model.getLensParameters()

# We'll store the density maps and the area that's considered in these
# variables.
area, densMaps = None, [ ]

try:
    # This is again a way to speed up the generation of the documentation:
    # if the density maps of the individual lenses have already been 
    # calculated, they can be retrieved from a file
    densMaps, area = pickle.loads(open("densmaps.bin", "rb").read())
except Exception as e:
    print(e)

    # If the density maps don't exist as a file yet, we'll calculate them.
    # For each lens in the CompositeLens, we'll request a calculation of
    # the density points in a certain region. This information is stored
    # so that the standard deviation can be calculated further down.
    for p in params:
        lensInfo = LI(p["lens"], size=120*ANGLE_ARCSEC)
        densMaps.append(lensInfo.getDensityPoints())
        area = lensInfo.getArea()

    # Save the calculated information to a file, to speed up the process
    # of (re-)generating the documentation.
    open("densmaps.bin", "wb").write(pickle.dumps((densMaps,area)))
    
# Calculate the standard deviation of all the individual density maps
# that were calculated
std = np.std(densMaps, 0)

# Create a DensInfo object that links the specified points, containing
# the standard deviation, to a plot area (dictionary with bottom-left and 
# top-right settings)
densInfo = DI(std, **area)

# Plot it in the same way we'd normally plot a GravitationalLens or LensInfo
# object.
plt.figure(figsize=(5,5))
plotutil.plotDensity(densInfo)
plt.gca().invert_xaxis()
plt.show()
