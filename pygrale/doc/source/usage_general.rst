.. _usage-module-general:

Inversion module 'general'
==========================

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

 - ``pointimageoverlap``: for strongly lensed point sources, this is a measure
   for how well the back-projected point images overlap in the source plane,
   as described in Zitrin et al (2010) [#2]_. By default the size of the
   area encompassing all backprojected points is used as a distance scale
   in determining the overlap, but the median of absolute deviations
   `MAD <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_ can be
   used as well (see the section about `'module parameters'`).

 - ``pointgroupoverlap``: an experimental fitness measure, which backprojects
   specific points and estimates the RMS in the image plane. Either the
   average source position is used, or each backprojected image is used
   separately as a source position estimate (see the section about `'module parameters'`).
   The magnification matrix is used to map differences in the source plane to 
   differences in the image plane.

 - ``extendedimagenull``: this fitness measure takes the null space into 
   account, the region where no images should be found, as described in
   Liesenborgs et al (2007) [#3]_. An extension was made so allow it to add
   a penalty when multiple images are generated when only a single image
   should be allowed.

 - ``pointimagenull``: similar to the null space fitness measure for extended
   images but now for point images, also described in Zitrin et al (2010) [#2]_.
   Here too the extension was made to allow singly imaged sources to be used
   as a contraint.

 - ``timedelay``: when time delay information is available for one or more
   multiply-imaged sources, a fitness measure will be calculated as described
   in Liesenborgs et al (2009) [#4]_ or as in Liesenborgs et al (2020) [#5]_.
   You can select which type to use using the `fitness_timedelay_type`
   parameters (see the section about `'module parameters'`).

 - ``causticpenalty``: also described in Liesenborgs et al (2009) [#4]_, it is
   also possible to add a fitness measure to try to prevent a caustic from
   intersecting a reconstructed source.

 - ``weaklensing``: if shear information is provided, a :math:`\chi^2`-like fitness 
   measure will be provided for this, as described in Liesenborgs et al (2020) [#5]_.
   By default, the shear data will be interpreted as averaged ellipticities, but
   if desired (for testing purposes for example), real shear, or actual reduced
   shear can be specified as well (see the section about `'module parameters'`).
   For each data point, a weight can be specified.

 - ``kappathreshold``: if it's certain that the convergence (:math:`\kappa`) values
   at certain locations should not exceed a certain threshold, this fitness
   measure is added.

 - ``bayesweaklensing``: TODO, bayesian, per galaxy

In general, when more than one fitness measure is used, at the end of the 
algorithm a set of mass maps is found in which no single one will be better
than another with respect to all fitness measures. While you can specify that
all of these solutions need to be retrieved, usually you'll just want
to keep a single solution. To choose this solution, the order of importance 
of the fitness measures can be left to the default or you can specify it 
yourself. See the section about `'module parameters'` for more information.

Input images data list
----------------------

In :func:`InversionWorkSpace.addImageDataToList <grale.inversion.InversionWorkSpace.addImageDataToList>`,
you specify the image type, as well as (optionally) extra parameters.
The images data set type name describes what kind of data
is contained in this entry, and can be one of the following:

 - ``pointimages``: the images data set describes point images that originate
   from a single source (also a point of course). Alternatively, extended
   images may be provided as well, but in this case point groups specifying
   which points are really images of the same source plane point must be
   present. In general in a strong
   lensing scenario there will be more than one image, but single images are
   allowed as well. If this is the case, the image can't be used directly 
   because it originates from a single point source by definition, but it can
   still be useful for a null space constraint to inform the algorithm that
   this source should not lie in a region that produces multiple images.

   For the overlap fitness, two options are possible: if the same ``groupname``
   is used for more multiple point images, the overlap of the points will be
   measured relative to the scale of only these back-projected points. By
   default, the scale of all backprojected points is used. If ``usescalefrom``
   is present, then for those points the scale set by another set of points,
   identified by a specific ``groupname``. 

   In case time delay information was added to the images data set, by 
   default a time delay fitness measure will be calculated as well. This can
   still be disabled by specifying a boolean extra parameter ``timedelay`` that
   is set to ``False``.

   The fitness measures to which this type of data is relevant, are 
   ``pointimageoverlap``, ``pointimagenull`` and ``timedelay``.

 - ``pointgroupimages``: an image data set with type is used in the
   ``pointgroupoverlap`` fitness measure, which calculates the RMS in the
   image plane. The input should either be point image data, or extended
   images in which point groups are used to indicate which points belong
   together.

   In case time delay information was added to the images data set, by 
   default a time delay fitness measure will be calculated as well. This can
   still be disabled by specifying a boolean extra parameter ``timedelay`` that
   is set to ``False``.

   The fitness measures to which this type of data is relevant, are
   ``pointgroupoverlap``, ``pointimagenull`` and ``timedelay``.

 - ``pointnullgrid``: the images data set should contain only one 'image', a
   triangulated grid, and is interpreted as a null space grid. It refers to 
   the previous images data set of ``pointimages`` type and the algorithm will 
   use this grid to estimate the number of extra images that are generated by 
   a trial solution. 
   
   One extra parameter can be specified: ``weight``, which should be a real and 
   positive value. The estimated amount of images is multiplied by this
   number, and can be used to adjust the relative importance of some null 
   spaces. This can be useful to provide extra weight to singly imaged sources
   because they may otherwise provide only a small penalty.

   The only fitness measure to which this type of data is relevant, is 
   ``pointimagenull``.

 - ``extendedimages``: in this case the images data set describes the extended 
   images originating from a single source. Therefore, in principle each image
   should consist of at least three points. A set of point images may be
   provided as well, but since each separate point image does not provide
   a distance scale to measure overlap with, the average size of the other
   backprojected (extended images) is used as the distance scale instead.
   A single image
   containing exactly one point is allowed as well. In this case it will not
   be used as a constraint regarding the overlap of back-projected images
   in the source plane, but it can be used as a constraint in the null space,
   indicating that the corresponding source plane position should produce only
   one image.
 
   For extended images, the extra parameters ``userectangles`` and 
   ``usepointgroups`` (both taking a boolean value) can be useful. By default,
   the overlap of surrounding back-projected rectangles is always used, but
   can be disabled with the first option. If point groups (points in different
   images that correspond to each other) are available in the images data set,
   they are used by default as well. The second option can disable their use.

   In case time delay information was added to the images data set, by 
   default a time delay fitness measure will be calculated as well. This can
   still be disabled by specifying a boolean extra parameter ``timedelay`` that
   is set to ``no``.

   The fitness measures to which this type of data is relevant, are 
   ``extendedimageoverlap``, `extendedimagenull```, ``pointimagenull`` and ``timedelay``.

 - ``extendednullgrid``: similar to ``pointnullgrid``, the images data set should 
   contain only one 'image', a triangulated grid, and is interpreted as a null 
   space grid. Here, typically the regions containing observed images or 
   regions that may harbor an as yet undetected image are removed from the 
   grid. The images data set refers to the previous images data set of 
   ``extendedimages`` type and the algorithm will use this grid to calculate the
   sizes of the regions in the image plane that contain additional images.
   
   One extra parameter can be specified: ``weight``, which should be a real and 
   positive value. The total size of the additional images is multiplied by 
   this number, and can be used to adjust the relative importance of some null
   spaces. This can be useful to provide extra weight to singly imaged sources
   because they may otherwise provide only a small penalty.

   The only fitness measure to which this type of data is relevant, is 
   ``extendedimagenull``.

 - ``sheardata``: use this type to provide shear measurements to the
   algorithm. The images data set should contain only one 'image', a set of 
   points in the image plane for which shear components have been
   specified. One extra parameter (a real number) called `threshold` must be
   provided: this contains a threshold value for :math:`|1-\kappa|`. Only when at
   a certain point the value for :math:`|1-\kappa|` exceeds the specified threshold,
   will it be included in the :math:`\chi^2`-like calculation.

   The only fitness measure to which this type of data is relevant, is 
   ``weaklensing``.

 - ``bayesellipticities``: TODO, galaxy ellipticities.
 
   The only fitness measure to which this type of data is relevant, is
   ``bayesweaklensing``.

 - ``kappathresholdpoints``: this data set should only contain a single 'image',
   a set of points at which the convergence :math:`\kappa` is calculated. The
   mandatory extra (real valued) parameter ``threshold`` specifies if a penalty
   is needed for a certain point: if its convergence is below the threshold,
   no penalty is added to the fitness measure, otherwise the amount by which
   the threshold is exceeded is added.

   The only fitness measure to which this type of data is relevant, is 
   ``kappathreshold``.

 - ``causticgrid``: when an images data set of this type is specified, it should
   contain a triangulated grid, which will be used to estimate the caustics.
   The algorithm will calculate the length of the caustics that intersect the 
   estimated source, and this source estimate is based on the previously
   encountered images data set marked as ``extendedimages``.
 
   The only fitness measure to which this type of data is relevant, is 
   ``causticpenalty``.

 - ``singlyimagedpoints``: the images data set of this type should contain only
   one 'image' entry, which may consist of several points. Each point is assumed
   to have only a single image, i.e. it should not lie in a multiply imaged
   region. These data are therefore intended to be used together with a
   null space grid of type `pointnullgrid`. The points could also be specified 
   separately as ``pointimages``, but this way is more convenient if several 
   points are needed for which the same null space grid can be used.

   The only fitness measure to which this type of data is relevant, is 
   ``pointimagenull``.


Module parameters
-----------------

Extra parameters for this module can be set using the `fitnessObjectParameters`
argument in e.g. :func:`inversion.invert <grale.inversion.InversionWorkSpace.invert>`.
The defaults can be obtained using the command
:func:`inversion.getDefaultModuleParameters <grale.inversion.getDefaultModuleParameters>`, 
and are listed in the following table:

================================================ =========================
**Parameter name**                               **Value**
------------------------------------------------ -------------------------
priority_causticpenalty                          100
priority_pointimagenull                          200
priority_extendedimagenull                       200
priority_pointgroupoverlap                       250
priority_pointimageoverlap                       300
priority_extendedimageoverlap                    300
priority_timedelay                               400
priority_weaklensing                             500
priority_bayesweaklensing                        500
priority_kappathreshold                          600
scalepriority_pointimageoverlap                  100
scalepriority_extendedimageoverlap               100
scalepriority_pointgroupoverlap                  200
scalepriority_pointimagenull                     -1
scalepriority_extendedimagenull                  -1
scalepriority_weaklensing                        300
scalepriority_timedelay                          -1
scalepriority_kappathreshold                     -1
scalepriority_causticpenalty                     -1
scalepriority_kappagradient                      -1
scalepriority_bayesweaklensing                   300
fitness_pointgroupoverlap_rmstype                'AllBetas'
fitness_pointimageoverlap_scaletype              'MinMax'
fitness_timedelay_type                           'NoSrc'
fitness_timedelay_nosrc_cutoff                   0.0
fitness_weaklensing_type                         'AveragedEllipticities'
fitness_bayesweaklensing_zdist_values            None
fitness_bayesweaklensing_zdist_range             None
fitness_bayesweaklensing_zdist_numsamples        16
fitness_bayesweaklensing_b_over_a_distribution   None
fitness_bayesweaklensing_sigmafactor             3.0
fitness_bayesweaklensing_sigmasteps              7
================================================ =========================

In case input is provided with ``pointgroupimages`` type, the ``pointgroupoverlap``
fitness calculation is used, which estimates the RMS in the image plane. It
does this by projecting the image points onto the source plane, determining
a source position based on these points, and using the magnification matrix
to convert differences in the source plane to difference in the image plane.
By default, each backprojected image point is used as a possible source
position, corresponding to the value ``AllBetas`` of ``fitness_pointgroupoverlap_rmstype``.
To use the averaged source position instead, you can set it to ``AverageBeta``.

If input is provided of ``pointimages`` type, the ``pointimageoverlap`` fitness
calculation will be activated. It projects the image points onto the source
plane, and uses the differences between the backprojected points to base the
fitness measure on. The distance scale with which these differences are measured
depends on all backprojected points and by default the size of the entire area
is used. This corresponds to the setting ``fitness_pointimageoverlap_scaletype``
to ``MinMax``. In case the `median of absolute deviations <https://en.wikipedia.org/wiki/Median_absolute_deviation>`_
should be used instead, this option can be set to ``MAD``.

The ``fitness_timedelay_type`` can also be ``Paper2009``, to use the older one from
the 2009 article [#4]_. TODO: describe ``"fitness_timedelay_nosrc_cutoff``

For the weak lensing fitness, the data is by default interpreted as (averaged)
ellipticity measurements. This corresponds the default value of ``AveragedEllipticities``
for ``fitness_weaklensing_type``. In case true shear is supplied, or the actual
reduced shear, the value can also be ``RealShear`` or ``RealReducedShear`` respectively.
This is mainly meant for testing purposes.

As the names suggest, the options that start with `priority_` describe priorities 
for fitness measures. These values do **not** have any effect on the way the 
genetic algorithm operates, only on the final solution that is chosen from the 
set of 'best' solutions.

Here, the lower the priority value, the more important it is considered to be,
so if for example both the ``extendedimageoverlap`` and ``extendedimagenull``
fitness measures are used based on the provided images data sets, there may
be more than one 'best' solution found by the genetic algorithm. You may
decide to save all these solutions by setting ``returnNds`` to ``True`` in
the :func:`invert <grale.inversion.InversionWorkSpace.invert>` function, but
typically you'll want to go on using a single solution. 
	
It is based on the priorities that were specified, that one solution will be 
chosen. Using the default settings in our example, first the solution(s) will 
be chosen with the lowest value for the ``extendedimagenull`` fitness value,
and then for ``extendedimageoverlap``. You can force this order to be turned
around by overriding some of these values in the ``fitnessObjectParameters``
argument.

As you can see, some priorities have the same value. In case two fitness
measures with the same priority are used, the number of images is typically
used as a tie breaker. For example, if you've specified both point images
and extended images as input, the fitness measures ``pointimageoverlap`` and
``extendedimageoverlap`` will be used. When choosing a single final solution,
the fitness value that corresponds to most images will have the best
priority for the default settings. So if there are more point images than
extended images, the point image criterion will be considered first.

.. rubric:: References

.. [#1] `A genetic algorithm for the non-parametric inversion of strong lensing systems <https://ui.adsabs.harvard.edu/abs/2006MNRAS.367.1209L/abstract>`_ 
.. [#2] `Full lensing analysis of Abell 1703: comparison of independent lens-modelling techniques <https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.1916Z/abstract>`_
.. [#3] `Non-parametric inversion of gravitational lensing systems with few images using a multi-objective genetic algorithm <https://ui.adsabs.harvard.edu/abs/2007MNRAS.380.1729L/abstract>`_
.. [#4] `Non-parametric strong lens inversion of SDSS J1004+4112 <https://ui.adsabs.harvard.edu/abs/2009MNRAS.397..341L/abstract>`_
.. [#5] `Extended lens reconstructions with grale: exploiting time-domain, substructural, and weak lensing information <https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.3253L/abstract>`_
