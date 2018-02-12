
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

**Status**

 * Creating various lenses, making plots etc should be fairly complete
 * Inversion code is starting to work from python as well, but needs examples
   and documentation
 * While several modules and classes are well documented, overall documentation
   and examples are still lacking

Contents
--------

.. toctree::
   :maxdepth: 2

   grale_constants.rst
   grale_cosmology.rst
   grale_feedback.rst
   grale_images.rst
   grale_lenses.rst
   grale_plotutil.rst
   grale_renderers.rst


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Installation
------------

The easiest way to install *pygrale* is using the `conda <http://conda.pydata.org/docs/intro.html>`_
tool of the `Anaconda Python <https://www.continuum.io>`_ distribution. Just run the command::

    conda install -c jori pygrale grale2modules

This should work on all platforms.

Overview
--------

Simulation
^^^^^^^^^^

TODO

Inversion
^^^^^^^^^

TODO

