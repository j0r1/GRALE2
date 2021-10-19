
GRALE2: pygrale
===============

GRALE is a project to allow you to simulate and invert
`gravitational lenses <https://en.wikipedia.org/wiki/Gravitational_lens>`_.
The `first generation <https://research.edm.uhasselt.be/jori/page/Physics/GraleV1.html>`_ of GRALE
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

Tutorial
--------

Continue to the :ref:`tutorial`.

.. _installation:

Installation
------------

There are two scripts available to help you set up everything. The easiest way to 
install *pygrale* is in a `conda <https://conda.io/en/latest/index.html>`_
environment of the `Anaconda Python <https://www.anaconda.com/>`_ distribution,
and the `BuildInCondaEnvironment.py <https://github.com/j0r1/GRALE2/tree/master/scripts>`_
script is intended to set this up automatically. This should work on Linux,
OS X and Windows.

The second one, `BuildAll.sh <https://github.com/j0r1/GRALE2/tree/master/scripts>`_,
is Linux specific (may work on OS X but is untested) and is mainly intended to 
get *pygrale* working on a supercomputer, but you can also use this script to setup a
`Google Colab <https://colab.research.google.com/>`_.

I'm very curious to see how well these scripts work, so feel free to
`contact me <mailto:jori.liesenborgs@gmail.com?subject=GRALE_Compilation_experience>`_ 
about your experiences with them (both failed and successful)!

For instructions on how to compile everything yourself, take a look at the
`GitHub <https://github.com/j0r1/GRALE2>`_ page.

Setting up a Colab
^^^^^^^^^^^^^^^^^^

Copy-paste and run this in a cell::

    !apt install libgsl-dev
    !apt install libopenmpi-dev
    !rm -f BuildAll.sh
    !wget https://raw.githubusercontent.com/j0r1/GRALE2/master/scripts/BuildAll.sh
    !NOQT=1 NOVENV=1 bash -e ./BuildAll.sh /usr/local/
    !/sbin/ldconfig

This takes a few minutes, but if everything was successful you can try out pygrale
in the Colab.

Setting up a conda environment with BuildInCondaEnvironment.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

I used to create ready-to-install packages, to allow one to run something like
`conda install pygrale`, but maintaining these for three different platforms
was not particularly interesting. This script is an alternative way to set up
everything in a conda environment: it will install some packages, and build the
rest from source.

The prerequisite is that you already have `Anaconda <https://www.anaconda.com/>`_
installed, as well as `git <https://git-scm.com/downloads>`_. For Windows, you 
also need to have the
`Visual Studio 2017 build tools <https://visualstudio.microsoft.com/thank-you-downloading-visual-studio/?sku=BuildTools&rel=15>`_
(see also `here <https://stackoverflow.com/questions/57795314/are-visual-studio-2017-build-tools-still-available-for-download>`_
for other links) installed.

The build procedure is slightly different for Windows and OS X/Linux.

Windows
,,,,,,,

To be able to compile the source code, you need to have the `Visual Studio 2017 
build tools` mentioned above installed. It has to be that particular version,
since we'll be using the `conda-forge <https://conda-forge.org/>`_ packages,
which use this version (see also `this <https://conda-forge.org/docs/maintainer/knowledge_base.html#local-testing>`_).

We'll need to have a command prompt with these VS 2017 tools as well as the conda toolset,
and to do this, go to your start menu and look for the `x64 Native Tools Command
Prompt for VS2017` (just start typing this). When you execute this, a command
line prompt will appear, and we need to activate Anaconda python in this. You'll
need to know in which location you have installed this, and then type something
like::

    c:\path\to\anaconda\Scripts\activate

After this, the ``conda`` command should be recognized. Then, use the ``cd`` command
to go to the path where the script `BuildInCondaEnvironment.py` is installed. This
can be either the path to the git repository if you've cloned that, or alternatively
the path in which you've downloaded that file from the `repository <https://github.com/j0r1/GRALE2/tree/master/scripts>`_::

    cd \path\to\script\folder

Just executing the script with ``python BuildInCondaEnvironment.py`` will show you
the instructions; it also mentions that you'll need to have the VS 2017 build tools
for example. As mentioned in the instructions, first create a conda environment in 
which everything GRALE related will be installed, something like::

    conda create -n mypygrale -c conda-forge -y python=3.8

As you can see here, we're specifying from the beginning that packages will be
used from the conda-forge channel. When this is finished, activate that environment::

    conda activate mypygrale

Then, decide on the directory where you want to store the downloaded source
code (git clones) as well as the intermediate compilation files, and execute::

    python BuildInCondaEnvironment.py c:\path\to\build\directory

If all goes well, this automatically clones the needed git repositories, compiles
the source code with the VS 2017 build tools, and installs everything in the
conda environment.

Once everything has been built, you don't need to run the `x64 Native Tools Command
Prompt for VS2017` anymore to activate the environment. You can just start an
Anaconda prompt from the start menu, and type ``conda activate mypygrale`` to
make everything available.

Linux & OS X
,,,,,,,,,,,,

After you've activated the Anaconda environment, so that the `conda` command
is available, first create a new environment in which all things for GRALE
will be installed::

    conda create -n mypygrale -c conda-forge -y python=3.8

As you can see here, we're specifying from the beginning that packages will be
used from the conda-forge channel. When this is finished, activate that environment::

    conda activate mypygrale

Use the ``cd`` command to go to the path where the script `BuildInCondaEnvironment.py`
is installed. This can be either the path to the git repository if you've 
cloned that, or alternatively the path in which you've downloaded that file
from the `repository <https://github.com/j0r1/GRALE2/tree/master/scripts>`_::

    cd /path/to/script/folder

Then, decide on the directory where you want to store the downloaded source
code (git clones) as well as the intermediate compilation files, and execute::

    python BuildInCondaEnvironment.py /path/to/build/directory

This will start with installing extra packages in the conda environment, and
you will notice that at a certain point it will ask you to deactivate and
reactivate the conda environment. This is to make sure that the compiler tools
that are installed as a conda package, will be set up correctly. After
reactivating the conda environment as instructed, just start the script again
in exactly the same way as before, and the process will continue where it
left off.

If all goes well, this automatically clones the needed git repositories, compiles
the source code, and installs everything in the conda environment.

Building with BuildAll.sh, e.g. on a supercomputer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `BuildAll.sh` script is intended for a Linux based environment, although
you might be able to use it on OS X as well. This does not make use of a conda
environment, but will create a standard Python
`virtual environment (virtualenv) <https://docs.python.org/3/tutorial/venv.html>`_,
in which several other packages will be installed using `pip <https://pypi.org/project/pip/>`_,
and in which the GRALE source code will be compiled and installed. To obtain
the script, either clone the git repository, or download it directly from the
`script directory <https://github.com/j0r1/GRALE2/tree/master/scripts>`_
(and set the executable flag with ``chmod u+x BuildAll.sh``)

Other than Python itself, the following things need to be available:

 - a C/C++ compiler
 - `CMake <https://cmake.org/>`_
 - `the git version control system <https://git-scm.com/>`_
 - `GSL, the GNU Scientific library <https://www.gnu.org/software/gsl/>`_

Not really required to get something to work, but definitely needed if you
want to use the compute resources from multiple nodes at the same time, is
an `MPI <https://en.wikipedia.org/wiki/Message_Passing_Interface>`_ implementation.

On a supercomputer, you'll typically need to activate these things using some
module system. For example, on one cluster I need to perform the following
commands to activate these ingredients::

     module load Python/3.6.4-intel-2018a    # A recent Python interpreter
     module load iimpi/2018a                 # The Intel compiler and Intel MPI tools
     module load CMake/3.13.3
     # git is available by default, no module needs to be loaded
     module load GSL/2.4-GCCcore-6.4.0

On another, the following modules need to be loaded::

    # A recent Python version is available by default
    module load gcc/8.2.0                    # The GNU C/C++ compilers
    module load cmake/3.10.2
    # git is available by default, no module needs to be loaded
    module load gsl/2.5
    module load ompi/4.0.0/gnu-8.2.0-centos7 # The OpenMPI MPI implementation


The `BuildAll.sh` script takes one argument, a directory which will be used to
set up the Python virtual environment in, as well as in which the source code
will be downloaded and the compiled code will be installed, e.g.::

    ./BuildAll.sh /path/to/new/virtualenv/directory

Some environment variables can be used to affect parts of the build process,
for a complete overview just run `BuildAll.sh` without any arguments. We'll
explore some of them next. On one system, I needed to make sure that the Intel compiler 
was used, and running

.. code-block::

    export CC=icc
    export CXX=icc
    
before executing the script made sure that this was the case. Similarly, on another
system I needed to enforce that `gcc` and `g++` were used as compilers, which could
be done by first setting

.. code-block::

    export CC=gcc
    export CXX=g++

Making the CMake build system detect the correct version of GSL sometimes also
required some extra work. In one case, the correct GSL library directory needed to
be specified, which can be done with the help of the `gsl-config` utility that comes 
with a GSL installation::

    export CMAKE_EXTRA_OPTS="-DADDITIONAL_LIBRARIES_CORE=-L`gsl-config --prefix`/lib"

On another system, all the paths to the GSL include files and libraries needed
to be specified explicitly::

    export CMAKE_EXTRA_OPTS="-DGSL_INCLUDE_DIR=/path/to/gsl/include -DGSL_CBLAS_LIBRARY=/path/to/libgslcblas.so -DGSL_LIBRARY=/path/to/libgsl.so"

By default, the build system will try to run as many things as possible in parallel,
to speed up the process. If you're running the script on a login node for example, that's
shared with other users, you may want to limit the number of parallel processes. This
can be done with, e.g.::

    export NUMCORES=4

The build script will try to detect if Qt5 is available, and if so, download the
compatible Qt interface for Python (`PyQt5 <https://riverbankcomputing.com/software/pyqt>`_).
This is needed to run the :ref:`GRALE Editor <graleeditor>`. Typically you don't
need this GUI program on the supercomputer however, just on your local machine.
You can bypass this step entirely, even if Qt5 is detected, by setting::

    export NOQT=1

When the script was able to successfully build everything, it ends with a line
that shows you how to enable this Python environment with the GRALE tools
installed, something like::

    source /path/to/new/virtualenv/directory/bin/activategrale

Note that before running this, e.g. when logging on to the supercomputer
again, or in a Slurm or TORQUE/PBS job description, you'll need to load
the same modules again (instead of the CMake one, that's only needed for
building the source code).

`Let me know <mailto:jori.liesenborgs@gmail.com?subject=Building_GRALE>`_
how it goes! I'd be happy to help out if you're having problems.

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
   grale_all.rst
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
   grale_multiplanecuda.rst
   grale_plotutil.rst
   grale_renderers.rst
   grale_util.rst
   grale_editor.rst
   inversion_module_usage.rst
   internals.rst
   debugging.rst
