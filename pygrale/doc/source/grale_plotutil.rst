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

.. autofunction:: plotDensityContours

.. autofunction:: plotDensityMixed3D2D

.. autofunction:: plotDensitiesAtImagePositions

.. autofunction:: plotImagePlane

.. autofunction:: plotAverageDensityProfile

.. autofunction:: plotIntegratedMassProfile

.. autofunction:: plotSubdivisionGrid

.. autofunction:: plotImagesData

.. autofunction:: plotShearComponents

Functions to create gnuplot plots
---------------------------------

.. autofunction:: plotDensityGnuplot

.. autofunction:: plotImagePlaneGnuplot

Functions to create FITS files
------------------------------

.. autofunction:: createEmptyFITS

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

.. autofunction:: getDensitiesAtImagePositions

.. autofunction:: mergeDensityMeasurementsAndAveragePositions

Classes to create animations
----------------------------

.. autoclass:: Animation
   :members:

   .. automethod:: __init__

.. autoclass:: NotebookAnimation
   :members:

