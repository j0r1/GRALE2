GRALE2: pygrale
===============

GRALE is a project to allow you to simulate and invert
[gravitational lenses](https://en.wikipedia.org/wiki/Gravitational_lens).
The [first generation](http://research.edm.uhasselt.be/jori/grale) of GRALE
consisted of C++ libraries and an interactive program called _GRALESHELL_
with which many useful commands could be executed without having to write
a program yourself using the C++ libraries.

The new generation _GRALE2_ consists of a trimmed down set of core libraries
that are still written in C++, together with [Python](http://www.python.org)
bindings to provide the same functionality: _pygrale_. The idea is that the core 
code is still written in the fast low-level C++ language, but to make everything
more easily useable and flexible, a Python interface is provided. The link
between the C++ code and Python could be provided thanks to the
[Cython](http://www.cython.org) project.

The basic idea is to put things that only depend on [GSL](https://www.gnu.org/software/gsl/),
[SerUt](https://github.com/j0r1/SerUt) and [ErrUt](https://github.com/j0r1/ErrUt) in the
core library, and the things that need [MOGAL](https://github.com/j0r1/MOGAL), and possibly 
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

### Using conda ###

The easiest way to install pygrale is using the [conda](http://conda.pydata.org/docs/intro.html)
tool of the [Anaconda Python](https://www.anaconda.com/) distribution. Just 
run the command:

    conda install -c jori pygrale grale2modules

This should work on all platforms.

### Manual build ###

Alternatively, you can build everything yourself. The first step is to compile
the libraries, for which you'll need a few dependencies:

 - Error utilities in [ErrUt](https://research.edm.uhasselt.be/jori/errut) 
   ([GitHub](https://github.com/j0r1/ErrUt) link)
 - Serialization utilities in [SerUt](https://research.edm.uhasselt.be/jori/serut)
   ([GitHub](https://github.com/j0r1/SerUt) link)
 - Platform independent network functions in [ENUt](https://research.edm.uhasselt.be/enut)
   ([GitHub](https://github.com/j0r1/ENUt) link)
 - The multi-objective genetic algorithm library [MOGAL](https://research.edm.uhasselt.be/jori/mogal)
   ([GitHub](https://github.com/j0r1/MOGAL) link)
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

If you're also interested in performing lens inversions, you'll need to compile
the inversion modules for the genetic algorithm as well. These are not compiled
automatically, but are provided in the subdirectory `inversion_modules`, where
the CMake build system can again be used to build them. You can either install
them on your system, or set the `GRALE2_MODULEPATH` environment variable to inform
the inversion scripts where to find them. This environment variable can also be
used in case they are not automatically detected after installation for example.



