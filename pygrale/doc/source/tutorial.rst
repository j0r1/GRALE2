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
methods that you can call for every one of the the available models. You can also
import a `LensTool <https://projets.lam.fr/projects/lenstool/wiki>`_ model using the
:func:`createLensFromLenstoolFile <grale.lenses.createLensFromLenstoolFile>` function.

The constructor of a specific lens model always takes two arguments. The first is the
angular diameter distance to the lens, the second describes the parameters for the
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
as well as the density map. In the example below, we use the 
:func:`plotDensityContours <grale.plotutil.plotDensityContours>` function
to create a contour map of the projected mass density of the lens.
The result is shown in the panel on the left.
The center and right panels of the same figure show the circularly averaged mass 
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
function. Note that actually a list of source shapes needs to be specified.
The example below illustrates this:

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

When creating a plot for a single source, or several sources at the same
distance/redshift, the last set angular diameter distances (or source 
redshift) is used. To plot multiple sources, the list of source shapes
needs to become a list of dictionaries containing the source shapes as
well as the distance information. The code below illustrates this.

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

Using a :class:`LensInfo <grale.plotutil.LensInfo>` object, you can also
simulate lensing by multiple gravitational lenses as different redshifts.
Instead of the :class:`LensPlane <grale.images.LensPlane>` and
:class:`ImagePlane <grale.images.ImagePlane>` classes, internally now
:class:`MultiLensPlane <grale.images.MultiLensPlane>` and
:class:`MultiImagePlane <grale.images.MultiImagePlane>` will be used.

For the `lens` parameter of the ``LensInfo`` object, you now need
to specify a list of tuples, where each of the tuples contains a
lens model and a redshift. To specify which cosmology must be used
to convert the redshifts to angular diameter distances, you can either
set the `cosmology` argument, or install a default cosmological model
using :func:`setDefaultCosmology <grale.cosmology.setDefaultCosmology>`
The example below illustrates this, more elaborate examples can be
found in the :ref:`notebooks` section.

.. literalinclude:: ex/example_plot_multiplane.py
   :language: python

.. plot:: ex/example_plot_multiplane.py

Regarding density plots, this cannot be done when the ``LensInfo``
object has been initialized for a multi-lensplane scenario. If you
want information about the mass distributions of the individual lenses,
you have to obtain the information for the individual lenses, and
combine them manually. This way, you can decide yourself whether
or not you need to take into account the lens effect of the foreground
lens, as is illustrated in the following example, which continues
from the SIS and SIE lenses in the previous one:

.. literalinclude:: ex/example_plot_multiplane_dens.py
   :language: python
   :lines: 39-

.. plot:: ex/example_plot_multiplane_dens.py

Modelling arbitrary mass distributions
--------------------------------------

There are several reasons you might want to be able to calculate the
lens effect for a mass distribution of which you don't know how to
expand it in terms of the basic :mod:`GravitationalLens <grale.lenses>`
derived classes: perhaps the lens shape is some analytic model and
either the lens effect is not known, or an appropriate model is
not available in Grale. Another reason might be that you want to know
the lens effect for some pixelized distribution, e.g. the light that's
observed in a cluster and stored in a FITS file.

For the pixelized case, in principle you could use the 
:class:`MultipleSquareLens <grale.lenses.MultipleSquareLens>` model
and add the effect for each individual pixel. This will certainly work,
and depending on the size of the pixels will create a good approximation
to the true, non-pixelized distribution. However, the number of pixels will become
very large rather quickly, causing the calculations for e.g. the
image plane (which itself is subdivided into a large number of locations)
to become very slow, even when you can accelerate the calculations
by using one of the :mod:`renderers <grale.renderers>`. Also, the final
model will lack some smoothness that may be desired, after all it is
built up from square pixels.

To reduce the complexity of the calculation and still obtain a good
model, it is possible to create a fit to the mass distribution using
Plummer basis functions of varying sizes. To do so, first you 
create a :mod:`grid <grale.grid>` that's based on the density you'd
like to model. Typically, this grid will be created so that the regions
that contain more mass will have a more narrow subdivision. Then,
using this grid as input, the :func:`fitMultiplePlummerLens <grale.grid.fitMultiplePlummerLens>`
function is called, which assigns Plummer basis functions to the
grid cells (with a width proportional to the cell size), fit the weights
of these basis functions to obtain a close match to the target density,
and finally returns the result as a :class:`MultiplePlummerLens <grale.lenses.MultiplePlummerLens>`
lens model.

The target density needs to be supplied to the ``fitMultiplePlummerLens``
function by supplying a callback function. Depending on what you'd like
to model, the :mod:`GridFunction <grale.gridfunction>` class can be
useful, for example to use a FITS file as input with the 
:func:`GridFunction.createFromFITS <grale.gridfunction.GridFunction.createFromFITS>`
function.

Examples of this fitting procedure can be found in the
`fittest.ipynb <_static/fittest.ipynb>`_ notebook.

How much the subdivision grid needs to be refined will depend on the
complexity of the mass distribution you're trying to model, and will
require some experimentation. In the various :func:`createSubdivisionGrid <grale.grid.createSubdivisionGrid>`-like
functions with which the grid will typically be created, there's an
argument called `keepLarger`, which defaults to ``False``. In my experience,
fitting is somewhat improved by setting this to ``True``; this causes
the original grid cell to be present in the subdivision grid even when
it's subdivided further.

Inversion
---------

While creating gravitational lens simulations is interesting in its own
right, usually one is far more interested in performing a gravitational lens
inversion. Below you can find some information to help get you started,
several `examples <https://github.com/j0r1/GRALE2/tree/master/inversion_examples>`_
are included in the source code archive to illustrate the procedure.

The origins of the inversion procedure start with 
`2006MNRAS.367.1209L <https://ui.adsabs.harvard.edu/abs/2006MNRAS.367.1209L/abstract>`_,
but over the years many additions were made, for example to include
so-called null-space information (where no images are present), or to
include time delay information. A recent overview of the entire procedure
can be found in 
`2020MNRAS.494.3253L <https://ui.adsabs.harvard.edu/abs/2020MNRAS.494.3253L/abstract>`_.

Overview
^^^^^^^^

You need to specify a few things to be able to start the inversion procedure:

 - the observations, usually multiple images of several sources, but could
   also include e.g. null space data, or weak lensing measurements;
 - a set of basis functions (type, location and initial shape/mass)
   for the lensing mass distribution, of which the weights need to be 
   determined to fit the observations specified as the input.

A genetic algorithm (GA) is then used as the optimization procedure, to
look for the weights of these basis functions. The GA looks *only* for
appropriate weights, so e.g. the position or orientation of a basis function
will not change.

To be able to model quite arbitrary mass distributions, one typically
starts with basis functions (the default type is a :class:`Plummer <grale.lenses.PlummerLens>`
model) laid out in a regular grid pattern. Then, based on a first run
of the GA, and therefore a first estimate of the lensing mass distribution,
one creates a grid in which cells are subdivided more finely in regions
where more mass is detected. This is illustated in the figure below:

.. plot:: ex/example_subdiv.py

Each cell will again correspond to a basis function of which the weight
needs to be optimized, and the GA is executed again to perform this optimization.
This subdivision procedure is typically repeated a number of times, with
more and more subdivisions/basisfunctions, until the added complexity
no longer provides an improvement in the retrieved solution.

Because no particular shape of the mass distribution is assumed, other
than the fact that it can be be built up using many simple basis functions,
the method is often termed 'free-form' or 'non-parametric'. You could
also place specific models, e.g. :class:`SIE <grale.lenses.SIELens>` models, at the
location of observed galaxies, and optimize their weights to fit the
observations. This is not unlike so-called parametric methods, but note that
the GA can only optimize the weights of the basis functions, it cannot
alter their location or rotation angle for example.

With an :class:`InversionWorkSpace <grale.inversion.InversionWorkSpace>`
object, it is quite straightforward to execute the needed steps:

 - first, you create an `InversionWorkSpace` where you specify the
   redshift of the lens, the size of region that will be used for the
   grid-based procedure, and the center of this region;

 - then, you :func:`add <grale.inversion.InversionWorkSpace.addImageDataToList>`
   contraints from observations to the workspace, specified as 
   :class:`ImagesData <grale.images.ImagesData>` objects;

 - you create a grid to lay out the basis functions, either a
   :func:`uniform <grale.inversion.InversionWorkSpace.setUniformGrid>`
   one, or a :func:`subdivision grid <grale.inversion.InversionWorkSpace.setSubdivisionGrid>`
   based on a previous lens model. Alternatively, you can also 
   manually :func:`add basis functions <grale.inversion.InversionWorkSpace.addBasisFunctions>`;

 - you :func:`run the GA <grale.inversion.InversionWorkSpace.invert>`
   to determine the weights of the basis functions, resulting in
   a mass model.

Based on the resulting mass model, you can create a new grid, start
the GA again, etc. This is the main procedure you'll find in the included
`examples <https://github.com/j0r1/GRALE2/tree/master/inversion_examples>`_.

Creating an InversionWorkSpace
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the procedure outlined above, during several refinement steps always the
same obvervational input data is used, and typically the same region is used
that is subdivided into a number of grid cells. To keep track of these
common settings, as well as to perform some other operations, an
:class:`InversionWorkSpace <grale.inversion.InversionWorkSpace>` instance
is created. 

For example, if a gravitational lens is located at a redshift of 0.4, 
the strong lensing mass should be retrieved in a 250 arcsec :sup:`2`
region, and a particular :class:`cosmological model <grale.cosmology.Cosmology>` 
is to be used, you'd do something like this::

   import grale.inversion as inversion
   import grale.cosmology as cosmology
   from grale.constants import *

   z_lens = 0.4
   cosm = cosmology.Cosmology(0.7, 0.27, 0, 0.73)
   iws = inversion.InversionWorkSpace(z_lens, 250*ANGLE_ARCSEC,
                                      cosmology=cosm)

Instead of specifying the cosmology in the constructor, you could also
first :func:`set a default <grale.cosmology.setDefaultCosmology>` one,
after which you can just omit the model in the construction of the
workspace, i.e. something like::

   cosm = cosmology.Cosmology(0.7, 0.27, 0, 0.73)
   cosmology.setDefaultCosmology(cosm)

   iws = inversion.InversionWorkSpace(z_lens, 250*ANGLE_ARCSEC)


Adding images
^^^^^^^^^^^^^

Next, we need to specify what the observational constraints are,
for example multiply imaged point sources. For such point image
data, it is common to have it listed in a text file, and there
exists a helper function :func:`readInputImagesFile <grale.images.readInputImagesFile>`
to process such a text file. It returns a list of dictionaries
with information about the redshift, the source identifier, and an
:class:`ImagesData <grale.images.ImagesData>` object (describing the
actual images), one for each multiply imaged source.
Depending on the format of such a text file, it is likely that you
need to specify how each line should be analyzed, the documentation
lists some examples. If no further processing is required, the
call would be something like this::

   import grale.images as images

   imgList = images.readInputImagesFile("inputpoints.txt", True)

For each image in the list, we then need to call the workspace's
:func:`addImageDataToList <grale.inversion.InversionWorkSpace.addImageDataToList>`.
In that call, we not only tell the workspace what the relevant redshift
is, but also what the type is of the data. In this case the data describes
point images, but it could also describe extended images or null space
information for example. In our example this would become::

   for i in imgList:
      iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")

For extended images it's also possible to read all the points from a
text file, but it's unlikely that such a file is available. Usually,
based on one or more FITS files of the sky region, you can use
the :ref:`GRALE Editor <tut-graleeditor>` tool to create these data sets,
and save them to a file. Again, there would be a different file for
each multiply imaged (extended) source. You could load and add these 
to the workspace like this::

   src1 = images.ImagesData.load("source1.imgdata")
   src2 = images.ImagesData.load("source2.imgdata")

   iws.addImageDataToList(src1, 2.5, "extendedimages")
   iws.addImageDataToList(src2, 1.5, "extendedimages")

A null space grid can be created using the helper function
:func:`createGridTriangles <grale.images.createGridTriangles>`,
or you could create it in the :ref:`GRALE Editor <tut-graleeditor>`. For
point images, the grid is typically a simple grid, but for extended
images the irrelevant regions should be cut out (the regions of
the images themselves, or a bright cluster galaxy for example).

To allow the inversion algorithm to figure out which null space data
belongs to which source, the data must be specified right after the
images themselves. For point images, as no holes are cut out there,
the same grid can typically be used for each source, and adding the
data could become something like this::

   bottomLeft = [ -200*ANGLE_ARCSEC, -200*ANGLE_ARCSEC]
   topRight = [ 200*ANGLE_ARCSEC, 200*ANGLE_ARCSEC]
   nullData = images.createGridTriangles(bottomLeft, topRight, 48, 48)

   for i in imgList:
      iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
      iws.addImageDataToList(nullData, i["z"], "pointnullgrid")

For extended images, as different regions *are* cut out for different
sources, a different file will need to be loaded. The code would then
become::

   src1 = images.ImagesData.load("source1.imgdata")
   null1 = images.ImagesData.load("null1.imgdata")
   src2 = images.ImagesData.load("source2.imgdata")
   null2 = images.ImagesData.load("null2.imgdata")

   iws.addImageDataToList(src1, 2.5, "extendedimages")
   iws.addImageDataToList(null1, 2.5, "extendednullgrid")
   iws.addImageDataToList(src2, 1.5, "extendedimages")
   iws.addImageDataToList(null2, 1.5, "extendednullgrid")

More information about the different types of input data that can
be specified this way, can be found in the `usage <./usage_general.html>`_
documentation.

**Whatever you do, before continuing, make sure that your input
data makes sense!** Create plots using
:func:`plotImagesData <grale.plotutil.plotImagesData>` or load
the data in the GRALE Editor.

Creating a grid
^^^^^^^^^^^^^^^
As explained before, the usual way to start a lens inversion is to
lay out basis functions in a regular grid pattern. In a next step,
lens plane basis functions (by default Plummer lenses) will be initialized to
have a width that's proportional to the size of a cell. To set up
such a regular grid, use :func:`setUniformGrid <grale.inversion.InversionWorkSpace.setUniformGrid>`,
for example::

   iws.setUniformGrid(15)

This would produce a regular 15x15 grid, with cells that cover the
region specified at the workspace's construction time. As the
documentation shows, by default some randomness will be applied to 
the grid center.

In case you already have a lens model, for example the result from
using this uniform grid (let's call it ``lens1``), you could 
create a new grid with cells that are finer in regions where 
there's more mass. This can be done with a call to
:func:`setSubdivisionGrid <grale.inversion.InversionWorkSpace.setSubdivisionGrid>`,
where you specify a range in which the number of cells should lie. 
As the 15x15 grid produces 225 cells, in a next step often the range 
from 300 to 400 is used, for example::

   iws.setSubdivisionGrid(lens1, 300, 400)

How many cells exactly this procedure would produce, depends on
the mass distribution. In a next iteration, with the optimization 
result from using this grid, you could for example create a new 
one as follows::

   iws.setSubdivisionGrid(lens2, 500, 600)

The grid that is created this way, can be obtained using 
:func:`InversionWorkSpace.getGrid <grale.inversion.InversionWorkSpace.getGrid>`,
and can be visualized easily using :func:`plotSubdivisionGrid<grale.plotutil.plotSubdivisionGrid>`.
This creates plots such as the ones shown above.

Running the inversion
^^^^^^^^^^^^^^^^^^^^^

The grid created in the previous step, determines the layout and widths
of the basis functions. A call to :func:`invert <grale.inversion.InversionWorkSpace.invert>`
determines the initial masses of these basis functions, and
subsequently starts the GA to figure out the weights of the basis
functions that are compatible with the observational constraints
provided to the workspace. In this `invert` call, a
so called 'inverter' can be specified, which allows you to speed
up the calculations using e.g. MPI. Alternatively, you can set a
default inverter using :func:`setDefaultInverter <grale.inversion.setDefaultInverter>`
and omit the parameter from the `invert` call. This would then
yield code like the following::

   inversion.setDefaultInverter("mpi")
   lens1, fitness, fitdesc = iws.invert(512)

The number 512 in the `invert` call, specifies the number of genomes/chromosomes/individuals
in the population of the GA. The function returns the :class:`GravitationalLens <grale.lenses.GravitationalLens>`
based lens model that is reconstructed, consisting of the basis functions 
with weights determined by the GA. This model can be 
:func:`saved <grale.lenses.GravitationalLens.save>` and later 
:func:`loaded <grale.lenses.GravitationalLens.load>` for further analysis.
Apart from this model, the `invert` call also returns the fitness value
that this model had in the GA, together with the type of the fitness measure
that was used (e.g `extendedimageoverlap` or `pointimageoverlap` depending
on the type of the provided input).

To assign initial masses to the basis functions, the mass of the
lens is estimated roughly from the provided images. If desired or needed,
this can be overridden using the `invert` call's ``massScale`` argument.
The default strategy is to simply divide this mass by the number of
grid cells, and assign each grid cell's basis function that share of the
mass. This will automatically place more mass in the more finely subdivided
regions. In case you feel this places too much mass in those regions, you
can set the ``rescaleBasisFunctions`` parameter to ``True``, causing these
weights to rescaled according to the size of the cell. In practice, the
default setting seems to provide better solutions faster.

Full control over which basis functions are placed where can be obtained
by using e.g. :func:`addBasisFunctions <grale.inversion.InversionWorkSpace.addBasisFunctions>`.
This allows you to use any lens model as a basis function. Instead of
`invert`, then the :func:`invertBasisFunctions <grale.inversion.InversionWorkSpace.invertBasisFunctions>`
call is required. As a side note, the `invert` call internally actually
calls this function once the basis functions have been derived from the
grid.

The ``massScale`` parameter also serves another purpose. Each genome in
the GA's population contains the weights of the basis functions for a
particular trial solution. As described in the articles, these weights
actually only determine the *shape* of the mass distribution. An overall 
scale factor will still be searched for, so that the resulting mass
distribution yields the best fitness value for that genome. This scale factor
will be probed in such a way that the total mass lies in a certain range 
around ``massScale``. The range can be controlled with the parameter
``massScaleSearchType``.

A special kind of basis function that can be enabled, is a mass sheet
basis function. It's special in the sense that it doesn't take part in
the scaling procedure mentioned above. To enable this, set the
``sheetSearch`` parameter to ``"genome"`` (the default is ``"nosheet"``).

To avoid having to specify the same arguments over and over in the `invert`
calls, you can use 
:func:`setDefaultInversionArguments <grale.inversion.InversionWorkSpace.setDefaultInversionArguments>`
to set the default values.

In case multiple types of constraints are used, e.g. multiple images as well
as null space information, a multi-objective GA will be used. Having multiple
constraints, and hence multiple fitness values, implies that there's in general
no longer a single best genome in each generation of the population, 
but a `non-dominated set`. 
The default behavior of the `invert` call is still to return a single
solution based on a priority of the individual fitness measures (can be 
controlled using the `priority` parameters, described in the
`usage <./usage_general.html>`_ documentation). To return the final non-dominated
set instead, set the ``returnNds`` parameter to ``True``. 

Processing the results
^^^^^^^^^^^^^^^^^^^^^^

The `invert` call returns a :class:`GravitationalLens <grale.lenses.GravitationalLens>`
based lens model, which you can then further analyze as described earlier in the
tutorial. For example, you can create plots of the mass distribution, or of the
critial lines and caustics. 

The workspace's function
:func:`backProject <grale.inversion.InversionWorkSpace.backProject>` can be
used to project the images back onto their respective source planes, 
using a particular lens model (doesn't even have to be one resulting from
the inversion routine). It returns
a list of `ImagesData` instances that have their point coordinates replaced
by the back-projected coordinates. Together with the
:func:`plotImagesData <grale.plotutil.plotImagesData>` function you could
create the plot below, showing the input images as well as the
back-projected ones for the first two lens models from the routine outlined
previously:

.. plot:: ex/example_backproj.py

The code to create such a plot looks like this:

.. literalinclude:: ex/example_backproj.py
   :language: python
   :lines: 15-34

With respect to re-lensing the source estimates, the functions
:func:`calculateRMS <grale.util.calculateRMS>` and
:func:`calculateImagePredictions <grale.util.calculateImagePredictions>`
may be of help.
These function reduce each observed image to a single point, which is
trivial for point images. In case you're dealing with extended images,
you can control how this needs to be done. 

In case you'd like to create a source shape from the backprojected images, the
:func:`createSourceFromImagesData <grale.images.createSourceFromImagesData>`
function may be of interest to you. The source shape this creates is just a simple 
filled polygon however, in case you'd like to see the actual lens effect using
one of the back-projected images as the source, you can calculate this with
the :ref:`GRALE Editor <tut-graleeditor>`.

You can always calculate the fitness value(s) a specific lens model (again,
doesn't need to be one that resulted from the inversion procedure) has
using :func:`calculateFitness <grale.inversion.InversionWorkSpace.calculateFitness>`.
For a lens model that was returned using the `invert` procedure, this should yield
roughly the same value as was returned by that function. Not precisely, as the
`calculateFitness` function uses code that's much different from the optimized
code in the inversion routine. The calculations can also take quite some time.

Note on reproducibility
^^^^^^^^^^^^^^^^^^^^^^^

In the procedure that's outlined above, there are actually two random number
generators at work, which is good to know if you'd like to be able to reproduce
your results exactly.

The first one is used in the creation of the grid: as
mentioned above, by default some randomness is added to the position of the
grid center. This randomness comes from Python's
`random module <https://docs.python.org/3/library/random.html>`_, so if
you want to make sure that the same grids are generated in different runs,
you can set the random number generator's seed using
`random.seed <https://docs.python.org/3/library/random.html#random.seed>`_.
Alternatively, you can save the random number generator's state using
`random.getstate <https://docs.python.org/3/library/random.html#random.getstate>`_
which you can then restore if needed using
`random.setstate <https://docs.python.org/3/library/random.html#random.setstate>`_.

Instead of working with Python's random number generator, you can also simply
save the generated grids. A grid can be obtained using the workspace's
:func:`getGrid <grale.inversion.InversionWorkSpace.getGrid>` function and can be
set again using :func:`setGrid <grale.inversion.InversionWorkSpace.setGrid>`.
To save or load very general Python objects, the easiest way is to
`pickle/unpickle <https://docs.python.org/3/library/pickle.html>`_ them.

To perform the inversion, so to run the GA, the Python code will actually
start a different program (written in C++), which has its own random number
generator. In the output of that program, you'll see a line like the following,
which allows you to find out which seed was used::

   GAMESSAGESTR:RNG SEED: 1581199847

To force a specific seed to be used when this inversion program is executed,
you can set the ``GRALE_DEBUG_SEED`` environment variable. For example,
to make sure that seed 12345 is used, one would run the code as follows::

   import os
   os.environ["GRALE_DEBUG_SEED"] = "12345"
   lens1, fitness, fitdesc = iws.invert(512)

Note that with this approach, the same environment variable will still be set
for other `invert` calls as well, which means that the same seed will be
used! To reproduce the results from other runs, you'd need to set
the seeds those different `invert` calls used, for example::

   os.environ["GRALE_DEBUG_SEED"] = "12345"
   lens1, fitness, fitdesc = iws.invert(512)

   iws.setSubdivisionGrid(lens1, 300, 400)
   
   os.environ["GRALE_DEBUG_SEED"] = "45678"
   lens2, fitness, fitdesc = iws.invert(512)

So, if you create an inversion script and want to be able to reproduce the
run exactly, you could print the random number generator's state at the
start, and perform a number of `invert` calls further on in the script,
as follows::

   import random
   print("RNG State:")
   print(random.getstate())

   ...

   iws = inversion.InversionWorkSpace( ... )

   ...

   iws.setUniformGrid(15)
   lens1, fitness, fitdesc = iws.invert(512)
   lens1.save("inv1.lensdata")

   iws.setSubdivisionGrid(lens1, 300, 400)
   lens2, fitness, fitdesc = iws.invert(512)
   lens2.save("inv2.lensdata")

   iws.setSubdivisionGrid(lens2, 500, 600)
   lens3, fitness, fitdesc = iws.invert(512)
   lens3.save("inv3.lensdata")

To be able to recreate this run exactly, you need to be able to inspect
all the output this script generated, so let's assume you've got that
in a file somehow. In that file, you can find the random number generator's
state (right below ``RNG State``, which we printed in the script) as well
as the seeds that were used in the external inversion program (those
``GAMESSAGESTR:RNG SEED:`` lines). Depending on that information, the
script to recreate the run would look something like this::

   import random
   random.setstate((3, (2869744549, 791801109, ... ), None))

   ...

   iws = inversion.InversionWorkSpace( ... )

   ...

   iws.setUniformGrid(15)
   os.environ["GRALE_DEBUG_SEED"] = "2551823485"
   lens1, fitness, fitdesc = iws.invert(512)
   lens1.save("inv1.lensdata")

   iws.setSubdivisionGrid(lens1, 300, 400)
   os.environ["GRALE_DEBUG_SEED"] = "1662455239"
   lens2, fitness, fitdesc = iws.invert(512)
   lens2.save("inv2.lensdata")

   iws.setSubdivisionGrid(lens2, 500, 600)
   os.environ["GRALE_DEBUG_SEED"] = "440472084"
   lens3, fitness, fitdesc = iws.invert(512)
   lens3.save("inv3.lensdata")


.. _tut-graleeditor:

GRALE Editor
------------

TODO: :ref:`GRALE Editor <graleeditor>`

Things to mention:

 - undo/redo
 - layers
 - draw image shapes, point groups, time delays
 - FITS background, align PNG overlay
 - null space grid, holes
 - export regions to file
 - backproject/re-lens
 - matching points/point groups, on backprojected image
