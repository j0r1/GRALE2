GRALE2: pygrale
===============

GRALE is a project to allow you to simulate and invert
[gravitational lenses](https://en.wikipedia.org/wiki/Gravitational_lens).
The [first generation](https://research.edm.uhasselt.be/jori/page/Physics/GraleV1.html) of GRALE
consisted of C++ libraries and an interactive program called _GRALESHELL_
with which many useful commands could be executed without having to write
a program yourself using the C++ libraries.

This new generation _GRALE2_ consists of a trimmed down set of core libraries
that are still written in C++, together with [Python](http://www.python.org)
bindings to provide the same functionality: _pygrale_. The idea is that the core 
code is still written in the fast low-level C++ language, but to make everything
more easily useable and flexible, a Python interface is provided. The link
between the C++ code and Python could be provided thanks to the
[Cython](http://www.cython.org) project.

The basic idea is to put things that only depend on [GSL](https://www.gnu.org/software/gsl/),
[SerUt](https://github.com/j0r1/SerUt) and [ErrUt](https://github.com/j0r1/ErrUt) in the
core library, and the things that need [EATk](https://github.com/j0r1/EATk), and possibly 
[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) in the inversion
library. The Python modules should only link against the core lib, thereby
avoiding a possible dependency on MPI. Interaction with e.g. MPI will be done 
using separate programs, that are communicated with using the [subprocess](https://docs.python.org/3/library/subprocess.html)
module.

Documentation
-------------

Documentation can be generated using [Sphinx](http://www.sphinx-doc.org/), or
can be viewed [online](http://research.edm.uhasselt.be/jori/grale2)

Installation
------------

The [documentation](http://research.edm.uhasselt.be/jori/grale2) describes how
you can use the provided [scripts](https://github.com/j0r1/GRALE2/tree/master/scripts)
to get your own copy working.

To really do everything yourself, the first step is to compile
the GRALE C++ libraries, for which you'll need a few dependencies:

 - Error utilities in [ErrUt](https://github.com/j0r1/ErrUt)
 - Serialization utilities in [SerUt](https://github.com/j0r1/SerUt)
 - The evolutionary algorithm toolkit [EATk](https://github.com/j0r1/EATk)
 - [GSL](https://www.gnu.org/software/gsl/), the GNU Scientific Library

Optionally, an MPI implementation can be very useful to speed up calculations; e.g.
with [OpenMPI](https://www.open-mpi.org/).

When these dependencies are satisfied, you can build the GRALE2 libraries: clone
the repository, and use [CMake](https://cmake.org/) to build and install them.

After these libraries are built _and_ installed, you can proceed to build the
Python bindings. To do so, you'll need [NumPy](http://www.numpy.org/) as well as
[Cython](http://cython.org). Just enter the `pygrale` directory and use the
`setup.py` script to build and install the Python modules, e.g. using

    python setup.py install

The following Python packages and programs will definitely be useful:

 - [SciPy](https://www.scipy.org/) is used for the integration needed to
   calculate angular diameter distances, as well as for some basic triangulation
   tasks.
 - [Astropy](http://www.astropy.org/) is used to read/write FITS files.
 - [Shapely](https://github.com/Toblerity/Shapely) is used to add a border
   to a polygon, when generating null space grids with holes in them for
   certain images.
 - When [PyQt5](https://www.riverbankcomputing.com/software/pyqt/download5) and
   [SIP](https://www.riverbankcomputing.com/software/sip/download) are available, 
   then the GRALE editor tool will be made available automatically. This package
   in turn depends on [Qt5](https://www.qt.io/).
 - Triangulations are made using the [triangle](http://www.cs.cmu.edu/~quake/triangle.html)
   program, which allows one to create constrained triangulations.

