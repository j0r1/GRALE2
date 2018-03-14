
GRALE2: pygrale
===============

GRALE is a project to allow you to simulate and invert
`gravitational lenses <https://en.wikipedia.org/wiki/Gravitational_lens>`_.
The `first generation <http://research.edm.uhasselt.be/jori/grale>`_ of GRALE
consisted of C++ libraries and an interactive program called *GRALESHELL*
with which many useful commands could be executed without having to write
a program yourself using the C++ libraries.

The new generation *GRALE2* consists of a trimmed down set of core libraries
that are still written in C++, together with `Python <http://www.python.org>`_
bindings to provide the same functionality: *pygrale*. The idea is that the core 
code is still written in the fast low-level C++ language, but to make everything
more easily useable and flexible, a Python interface is provided. The link
between the C++ code and Python could be provided thanks to the
`Cython <http://www.cython.org>`_ project.

This documentation is about the Python bindings only.

.. _installation:

Installation
------------

The easiest way to install *pygrale* is using the `conda <http://conda.pydata.org/docs/intro.html>`_
tool of the `Anaconda Python <https://www.continuum.io>`_ distribution. Just run the command::

    conda install -c jori pygrale grale2modules

in your environment. This should work on all platforms.

For instructions on how to compile everything yourself, take a look at the
`GitHub <https://github.com/j0r1/GRALE2>`_ page.

Tutorial
--------

Continue to the :ref:`tutorial`.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 1
   :hidden:
   
   tutorial.rst
   example_notebooks.rst
   grale_constants.rst
   grale_contourfinder.rst
   grale_cosmology.rst
   grale_feedback.rst
   grale_grid.rst
   grale_gridfunction.rst
   grale_images.rst
   grale_inversion.rst
   grale_inverters.rst
   grale_lenses.rst
   grale_multiplane.rst
   grale_plotutil.rst
   grale_renderers.rst


