"""This is a helper module that imports many of the grale modules,
sets various :mod:`renderers <grale.renderers>` to the 'threads' 
implementation, sets the 
:func:`default plotting angular scale <grale.plotutil.setDefaultAngularUnit>`
to arc seconds, and defines a function `V(x,y)` that creates a numpy 
array with two entries.
"""
from . import lenses
from . import renderers
from . import inversion
from . import inverters
from . import feedback
from . import images
from . import plotutil
from . import cosmology
from . import multiplane
from . import util
from . import grid
from grale.constants import *
import matplotlib.pyplot as plt
import numpy as np

V = lambda x, y: np.array([x, y], dtype=np.double)
renderers.setDefaultLensPlaneRenderer("threads")
renderers.setDefaultMassRenderer("threads")
print("Set 'threads' as default renderer for lensplane and mass density")

plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)
print("Set default angular unit in plotting to arcsec")
