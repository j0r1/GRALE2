.. _usage-module-general:

Inversion module 'general'
==========================

TODO: migrate usage info here

This lens inversion module can handle a variety of input data: image data,
null space grids, time delay info, etc. Based on the precise input that was
provided, one or more fitness measures will be used and the genetic
algorithm will try to optimize these. The currently implemented fitness 
measures are:

 - ``extendedimageoverlap``: for multiply imaged extended sources, this fitness
   measure calculates the fractional overlap of back-projected images as 
   described in Liesenborgs et al (2006) [#1]_. Point sources can be included
   as well, in that case the average size of the extended images is used
   as the distance scale for determining how well the back-projected points
   overlap.



UNDER CONSTRUCTION

.. rubric:: References

.. [#1] `A genetic algorithm for the non-parametric inversion of strong lensing systems <https://ui.adsabs.harvard.edu/abs/2006MNRAS.367.1209L/abstract>`_ 
.. [#2] `Full lensing analysis of Abell 1703: comparison of independent lens-modelling techniques <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.1916Z/abstract>`_
.. [#3] `Non-parametric inversion of gravitational lensing systems with few images using a multi-objective genetic algorithm <https://ui.adsabs.harvard.edu/abs/2007MNRAS.380.1729L/abstract>`_
.. [#4] `Non-parametric strong lens inversion of SDSS J1004+4112 <https://ui.adsabs.harvard.edu/abs/2009MNRAS.397..341L/abstract>`_
