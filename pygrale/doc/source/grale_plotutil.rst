grale.plotutil
==============

.. automodule:: grale.plotutil

Defaults
--------
.. autofunction:: setDefaultAngularUnit

.. autofunction:: getDefaultAngularUnit

DensInfo and LensInfo
---------------------

.. autoexception:: LensInfoException

.. autoclass:: DensInfo
   :members:

   .. automethod:: __init__

.. autoclass:: LensInfo
   :members:

   .. automethod:: __init__

PlotException
-------------
.. autoexception:: PlotException

Functions to create matplotlib plots
------------------------------------

.. autofunction:: plotDensity

.. autofunction:: plotImagePlane

.. autofunction:: plotAverageDensityProfile

.. autofunction:: plotIntegratedMassProfile

.. autofunction:: plotSubdivisionGrid

.. autofunction:: plotImagesData

Functions to create gnuplot plots
---------------------------------

.. autofunction:: plotDensityGnuplot

.. autofunction:: plotImagePlaneGnuplot

Functions to create FITS files
------------------------------

.. autofunction:: plotDensityFITS

.. autofunction:: plotImagePlaneFITS

.. autofunction:: calculateDeflectionAndDerivativesForFITS

Interactive plotting
--------------------

.. autofunction:: plot3DInteractive

.. autofunction:: plotDensityInteractive

Helper functions
----------------

.. autofunction:: estimatePlotScale

.. autofunction:: quickLensInfo

Classes to create animations
----------------------------

.. autoclass:: Animation
   :members:

   .. automethod:: __init__

.. autoclass:: NotebookAnimation
   :members:

