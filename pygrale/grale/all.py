from . import lenses
from . import renderers
from . import inversion
from . import feedback
from . import images
from . import plotutil
from . import cosmology
from . import multiplane
from grale.constants import *
import matplotlib.pyplot as plt
import numpy as np

V = lambda x, y: np.array([x, y], dtype=np.double)
renderers.setDefaultLensPlaneRenderer("threads")
renderers.setDefaultMassRenderer("threads")
print("Set 'threads' as default renderer for lensplane and mass density")

plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)
print("Set default angular unit in plotting to arcsec")
