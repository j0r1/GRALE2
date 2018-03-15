.. _tutorial:

Tutorial
========

After following the :ref:`installation` instructions, you'll obviously be
eager to get started with modeling and inversion of gravitational lens systems.
Instead of having to go through all available modules, classes and functions
first, this tutorial aims to familiarize you with the core ideas and procedures.
You can also take a look at various :ref:`notebooks` to see them in
action immediately.

Modeling gravitational lenses
-----------------------------

One of the core parts of the software is the :class:`GravitationalLens <grale.lenses.GravitationalLens>`
class. You won't create an instance of this class directly, that's prohibited, instead
you'll create a lens model using one of the derived classes, like :class:`SISLens <grale.lenses.SISLens>`
for lensing by a Singular Isothermal Sphere, or :class:`NFWLens <grale.lenses.NFWLens>` for
lensing by a mass distribution with a spherical Navarro-Frenk-White profile. 
All available models are shown in the documentation of the :mod:`grale.lenses` module;
the :class:`GravitationalLens <grale.lenses.GravitationalLens>` class contains the
methods that you can call for every one of the the available models.

The constructor of a specific lens model always takes two arguments. The first is the
angular diameter distance to the lens, the second describes the parameter for the
lens model. To convert the redshift of a lens to an angular diameter distance,
you can use an instance of the :class:`Cosmology <grale.cosmology.Cosmology>`
class, which destribes a particular cosmological model. Most often, various
quantities are expressed in SI units. In :mod:`grale.constants` there are various
constants that help you express values in these units, or to convert values to
other units.

The following example illustrates the creation of a lens model:

.. literalinclude:: ex/example_lens.py
   :language: python
   :lines: 1-17

In this example, we created a single SIS lens, which is centered at the origin.
There are many other available models, which are also typically centered at the
origin. To make it possible to use a different center, and to even combine multiple
lens models into a single, more complex model, we can use a :class:`CompositeLens <grale.lenses.CompositeLens>`.
Continuing from the previous example, we could change the center of the SIS lens,
and add a Singular Isothermal Ellipsoid at another location, even rotating it
over a specific angle (as the :class:`documentation <grale.lenses.CompositeLens>`
mentions, this angle is interpreted as degrees).

.. literalinclude:: ex/example_lens.py
   :language: python
   :lines: 19-28

You can save save a lens model to a file using the aptly named :func:`save <grale.lenses.GravitationalLens.save>`
function, and load the file again using the :func:`load <grale.lenses.GravitationalLens.load>` 
function. The binary data that is loaded or saved, can also be used directly through
the :func:`fromBytes <grale.lenses.GravitationalLens.fromBytes>` and
:func:`toBytes <grale.lenses.GravitationalLens.toBytes>` functions. This can help
to bypass having to write a lens model to a file first, e.g. when downloading
a model from the web, as the following example illustrates by loading the model
obtained by inverting CL0024:

.. literalinclude:: ex/example_lens_dl.py
   :language: python


Plotting lens systems
---------------------

The :mod:`plotutil <grale.plotutil>` module contains various functions to create plots, 
using `matplotlib <https://matplotlib.org/>`_, `gnuplot <http://gnuplot.info/>`_
or even an interactive 3D plot when using a `Jupyter <http://jupyter.org/>`_
notebook.

Renderers and feedback
^^^^^^^^^^^^^^^^^^^^^^

By default, calculations will be done on a single CPU core. For a complex lens
model, this can become quite slow, and various so called :mod:`renderers <grale.renderers>`
are available to speed up this process. There's one that can use the different
cores on your computer, and even one that can use different computers when the
`MPI <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_ system (e.g.
on a supercomputer) is available. For some calculations, `OpenCL <https://www.khronos.org/opencl/>`_
can use your GPU to speed up calculations, but this is not available for all
lens types. To get some feedback during such longer calculations, the
:mod:`feedback <grale.feedback>` module can be useful.

Many of the functions in the :mod:`plotutil <grale.plotutil>` module accept
`renderer` and `feedbackObject` arguments. While you can specify them in every
function used, it will probably be less cumbersome to simply set the default
to be used, using :func:`setDefaultMassRenderer <grale.renderers.setDefaultMassRenderer>`,
:func:`setDefaultLensPlaneRenderer <grale.renderers.setDefaultLensPlaneRenderer>`
and :func:`setDefaultFeedback <grale.feedback.setDefaultFeedback>`, for 
example:

.. literalinclude:: ex/example_rendererdefaults.py
   :language: python


Single lens plane
^^^^^^^^^^^^^^^^^

In the case of a single lens plane, we can create a plot of the image
plane and projected mass density with e.g. the :func:`plotImagePlane <grale.plotutil.plotImagePlane>`
:func:`plotDensity <grale.plotutil.plotDensity>` functions. The lens model
can be passed as the first argument, which will cause a relevant area to
be :func:`estimated <grale.plotutil.estimatePlotScale>`, which in turn will
be used for the plots. In the image plane case, a :math:`D_{ds}/D_s` fraction
of 1 will be assumed.

To avoid using radians as the unit on the plot axes, we can specify a particular
unit to be used with the `angularUnit` argument of the functions.
Alternatively, instead of specifying this argument with every single function, we
can also set the default unit using :func:`setDefaultAngularUnit <grale.plotutil.setDefaultAngularUnit>`
function. Note that this only converts values when plotting, it does not change
the way function parameters are interpreted for example.

The code below illustrates this, for the ``combLens`` model that was created in
a previous example:

.. literalinclude:: ex/example_plot.py
   :language: python
   :lines: 2-

This produces a matplotlib plot containing the two panels shown below, of which the left one
shows the image plane (for :math:`D_{ds}/D_s = 1`) and the right one shows the
density. In the left panel, a red line is used for the critical line and a blue
one for the caustic. In the right panel,
the `vmax` argument of the underlying `imshow <https://matplotlib.org/api/_as_gen/matplotlib.pyplot.imshow.html>`_
function is used to cut off the color scale at :math:`10 kg/m^2` to better show
the combination of the circularly symmetric SIS lens and the elliptic SIE lens.

.. plot:: ex/example_plot.py

To better control what to plot, e.g. the specific area or source distances,
instead of using a lens model as the first argument we can create a 
:class:`LensInfo <grale.plotutil.LensInfo>` instance for this situation.
Apart from more control about what to plot, this instance will also be used
to cache results, so that creating a similar plot may not require everything
to be recalculated. This is especially useful as the lens model becomes more
complex. When specifying a lens instead of a :class:`LensInfo <grale.plotutil.LensInfo>`
object, internally a ``LensInfo`` instance is still generated, and in fact, this
instance is typically returned by the functions. So in our previous example
we could have obtained this object in the following way:

.. code-block:: python

    lensInfo = plotutil.plotImagePlane(combLens)

The following example provides a more extensive example for plotting image
plane and mass density information, for the CL0024 lens from an earlier example.
In this example, we also set the renderers to ``threads``, so that all available
cores on a computer would be used. The plots that are generated by this example
are shown below the code (which continues where the CL0024 example left off):

.. literalinclude:: ex/example_plot_cl0024.py
   :language: python
   :lines: 2-

.. plot:: ex/example_plot_cl0024.py

The :class:`LensInfo <grale.plotutil.LensInfo>` object stores various 
properties, like the :class:`ImagePlane <grale.images.ImagePlane>` and
:class:`LensPlane <grale.images.LensPlane>` instances that were calculated,
as well as the density map. In the example below, we use this stored information
to generate a contour map using matplotlib's `contour <https://matplotlib.org/api/_as_gen/matplotlib.pyplot.contour.html>`_
function. Note that in this case we're using the raw matplotlib functions, and
our own plot scale factor will not be taken into account automatically.
For this reason, we convert to arcsec units manually. The result is shown in the
panel on the left.
The center and right panels of the figure above show the circularly averaged mass 
density profile and integrated mass profile respectively, using the 
:func:`plotAverageDensityProfile <grale.plotutil.plotAverageDensityProfile>`
and :func:`plotIntegratedMassProfile <grale.plotutil.plotIntegratedMassProfile>`
functions. 

.. literalinclude:: ex/example_plot_cl0024_b.py
   :language: python
   :lines: 10-

.. plot:: ex/example_plot_cl0024_b.py

Being able to plot the critical lines and caustics for a lensing scenario is
of course nice, but typically we're also interested in an image configuration
that corresponds to a specific source position. For this, we can create a
source shape using a class derived from :class:`SourceImage <grale.images.SourceImage>`,
and pass this to e.g. the :func:`plotImagePlane <grale.plotutil.plotImagePlane>`
function. Note that actually a list of source shapes needs to be specified; this
way you could also specify multiple sources (for this simple usage, they need to
be at the same redshift though). The example below illustrates this:

.. literalinclude:: ex/example_plot_cl0024_c.py
   :language: python
   :lines: 10-

.. plot:: ex/example_plot_cl0024_c.py

You could even create an animation to see what the effect is of different
source positions, using :class:`Animation <grale.plotutil.Animation>` or
:class:`NotebookAnimation <grale.plotutil.NotebookAnimation>`. 
Because the deflection calculations are cached, calculating the lens effect
at different source positions or even different source redshifts can be
done relatively quickly. The example
below is created in the notebook `fittest.ipynb <_static/fittest.ipynb>`_,
and shows the lens effect of a very strange mass distribution on a circular
source.

.. raw:: html

   <video src="_static/gralesource.mp4" width="50%" controls></video>

Creating plots for multiple sources at different redshifts takes a bit more
work. The most control is obtained by getting the :class:`ImagePlane <grale.images.ImagePlane>` 
instance using :func:`LensInfo.getImagePlane <grale.plotutil.LensInfo.getImagePlane>`,
from which you can obtain the :func:`caustics <grale.images.ImagePlane.getCaustics>`
and :func:`critical lines <grale.images.ImagePlane.getCriticalLines>`, and you
can render :func:`sources <grale.images.ImagePlane.renderSources>` and
:func:`images <grale.images.ImagePlane.renderImages>` which you can combine
with other sources or images. 

If a matplotlib figure is what you're interested in, the :func:`plotImagePlane <grale.plotutil.plotImagePlane>`
function also allows for some flexibility to combine the results for
different sources. The key is to make use of the `processRenderPixels`
argument, which should then be a function that will receive the NumPy
array that's about to be rendered as its argument. The return value of
this function will be what's actually plotted. The code below illustrates
this approach:

.. literalinclude:: ex/example_plot_cl0024_d.py
   :language: python
   :lines: 22-

.. plot:: ex/example_plot_cl0024_d.py

If you don't just want to plot the mass density of a gravitational lens,
but you want to look at the difference between lenses, or even the standard
deviation of a whole set of lenses, the :class:`DensInfo <grale.plotutil.DensInfo>`
class can come in handy. It accepts a 2D NumPy array of values, where each value is
considered to be the mass density at a specific point. These points are arranged regularly
between a bottom-left and top-right corner that also need to be specified.
For a specific lens, this kind of information can be obtained using a call to
:func:`LensInfo.getDensityPoints <grale.plotutil.LensInfo.getDensityPoints>`.

In the example below, we use the CL0024 model again. This is actualy an average
of several different individual lens models, and apart from the average mass
density associated to the lens, we're typically also interested in the way these
individual models vary. The standard deviation of all these maps can give us
some insight into this, and it's precisely this what is calculated in the
example:

.. literalinclude:: ex/example_plot_stddev.py
   :language: python
   :lines: 2-

.. plot:: ex/example_plot_stddev.py

Multiple lens planes
^^^^^^^^^^^^^^^^^^^^

:class:`grale.images.MultiImagePlane`,
:class:`grale.images.MultiLensPlane`




Fitting to a mass distribution
------------------------------

TODO: :mod:`grale.grid`, :mod:`grale.gridfunction`

Inversion
---------

TODO: :class:`grale.inversion.InversionWorkSpace`, :mod:`grale.inverters`,
:class:`images.ImagesData`, :func:`grale.images.readInputImagesFile`

Overview
^^^^^^^^

TODO

Adding images
^^^^^^^^^^^^^

TODO

Creating a grid
^^^^^^^^^^^^^^^

TODO

Running the inversion
^^^^^^^^^^^^^^^^^^^^^

TODO

Processing the results
^^^^^^^^^^^^^^^^^^^^^^

TODO

Examples
^^^^^^^^

TODO


