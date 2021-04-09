
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
get *pygrale* working on a supercomputer.

I'm very curious to see how well these scripts work, so feel free to
`contact me <mailto:jori.liesenborgs@gmail.com?subject=GRALE_Compilation_experience>`_ 
about your experiences with them (both failed and successful)!

For instructions on how to compile everything yourself, take a look at the
`GitHub <https://github.com/j0r1/GRALE2>`_ page.

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

Building on a supercomputer with BuildAll.sh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TODO

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
