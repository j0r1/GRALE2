grale.images
============

.. automodule:: grale.images

Functions
---------

.. autofunction:: centerOnPosition

.. autofunction:: uncenterOnPosition

.. autofunction:: readInputImagesFile

.. autofunction:: getDefaultLineAnalyzer

.. autofunction:: setDefaultLineAnalyzer

.. autofunction:: hoursMinutesSecondsToDegrees

.. autofunction:: degreesMinutesSecondsToDegrees

.. autofunction:: enlargePolygon

ImagesData
----------

.. autoclass:: ImagesData
   :members:

   .. automethod:: __init__

.. autoexception:: ImagesDataException

.. autofunction:: createGridTriangles

LensPlane & ImagePlane
----------------------

.. autoclass:: LensPlane
   :members:

   .. automethod:: __init__

.. autoexception:: LensPlaneException

.. autoclass:: ImagePlane
   :members:

   .. automethod:: __init__

.. autoexception:: ImagePlaneException

.. _sourceshapes:

Source shapes
-------------

.. autoclass:: SourceImage
   :members:

.. autoclass:: CircularSource
   :members:

   .. automethod:: __init__

.. autoclass:: EllipticalSource
   :members:

   .. automethod:: __init__

.. autoclass:: PolygonSource
   :members:

   .. automethod:: __init__

.. autoclass:: DiscreteSource
   :members:

   .. automethod:: __init__

.. autoclass:: PointSource
   :members:

   .. automethod:: __init__

.. autoexception:: SourceImageException

