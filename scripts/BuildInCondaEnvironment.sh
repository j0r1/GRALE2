#!/bin/bash -e

usage() {
cat << EOF

This script attempts to compile PyGrale in an existing conda environment,
where the specified directory stores the source code and build files.

Start with creating a new conda environment, and activate it:

    conda create -y mypygrale
    conda activate mypygrale

Then, running the script should install/build everything in this environment.

Of course: use at your own risk

Environment variables that are used:

 - NUMCORES: if set, 'make' will be started to use this amount of
   cores for the build (with '-j' option). If not set, the number
   of cores will be detected using the python multiprocessing module

 - NOENVCHECK: don't check the conda environment

 - NOPULL: skip git pull to update used repositories

 - NOCONDAINSTALL: assume all needed dependencies are available,
   don't run 'conda install' to add any packages.

 - UNAME: this is detected automatically by running the 'uname'
   command, but this environment variable overrides it. This is used
   to detect the correct compiler package for your platform.

 - PYVER: Python version to install in this environment, e.g. '3.7'

 - NOGSL: if you don't want to install a GSL package using anaconda,
   e.g. if you want to use a different one, you can set this.

 - NOMPI: if you don't want to install the openmpi package using
   anaconda, you can set this.

 - CMAKE_EXTRA_OPTS: extra options that are passed to each 'cmake'
   command. I use it to set 'ADDITIONAL_LIBRARIES_CORE' for GRALE2
   to be able to use the correct GSL libraries, with something like

    CMAKE_EXTRA_OPTS=-DADDITIONAL_LIBRARIES_CORE=-L\`gsl-config --prefix\`/lib

   (a newer GSL library needs to be used, enabled using a 'module load'
   command, but the system-wide one is used in the linking process. Seems
   similar as in e.g. https://cmake.org/pipermail/cmake/2011-June/044790.html,
   where the library's path is removed from the link command because it
   is also set in \$LIBRARY_PATH)

EOF
}

if [ "$#" != 1 ] ; then
	echo "Please specify a directory to store the build files"
	usage
	exit -1
fi

if ! git --version >/dev/null 2>&1 ; then
	echo "Need git"
	exit -1
fi

if ! conda --version >/dev/null 2>&1 ; then
	echo "Need active conda environment - 'conda' command is not available"
	exit -1
fi

if [ -z "$CONDA_PREFIX" ] ; then
	echo "CONDA_PREFIX environment variable is not set, is conda environment active?"
	exit -1
fi

if [ -z "$NOENVCHECK" ] ; then
	if [ -z "$CONDA_DEFAULT_ENV" ] ; then
		echo "No conda default environment could be determined"
		exit -1
	fi

	if [ "$CONDA_DEFAULT_ENV" = "base" ] ; then
		echo "Refusing to install in base conda environment - set NOENVCHECK to override"
		exit -1
	fi
	echo
	echo "Will install in conda environment '$CONDA_DEFAULT_ENV'"
	echo
	sleep 1
fi

STARTDIR=`pwd`
D="$1"
if ! [ -e "$D" ] ; then
	mkdir "$D"
fi
cd "$D"
D=`pwd`

if [ -z "$NUMCORES" ] ; then
	NUMCORES=`python -c "import multiprocessing;print(multiprocessing.cpu_count())"`
fi

if [ -z "$UNAME" ] ; then
	UNAME=`uname`
fi

if ! [ -z "$PYVER" ] ; then
	PYVER="=$PYVER"
fi

if [ -z "$NOCONDAINSTALL" ] ; then
	
	OWNPACKS="errut serut enut mogal triangle"
	if [ -z "NOGSL" ] ; then
		OWNPACKS="$OWNPACKS gsl"
	fi
	if [ -z "NOMPI" ] ; then
		OWNPACKS="$OWNPACKS openmpi"
	fi

	conda install -y -c jori $OWNPACKS

	CONDAPACKS="python$PYVER ipython jupyter astropy sip pyqt cython numpy scipy matplotlib shapely PyOpenGL ipywidgets cmake"
	if [ "$UNAME" = "Darwin" ] ; then
		CONDAPACKS="$CONDAPACKS clangxx_osx-64"
	elif [ "$UNAME" = "Linux" ] ; then
		CONDAPACKS="$CONDAPACKS gxx_linux-64"
	else
		echo "Unexpected uname of $UNAME, expecting Darwin or Linux. You can set the UNAME environment variable to force a choice"
		exit -1
	fi

	echo "Installing $CONDAPACKS"
	conda install -y $CONDAPACKS

fi

for p in GRALE2 ; do
	cd "$D"
	if ! [ -e $p ] ; then
		echo "Cloning $p"
		git clone https://github.com/j0r1/$p
		cd $p
	else
		echo "Updating $p"
		cd $p
		if [ -z "$NOPULL" ] ; then
			git pull
		else
			echo "Environment variable NOPULL was set, skipping git pull"
		fi
	fi

	echo "Building $p"
	if ! [ -e "build" ] ; then
		mkdir build
	fi
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" -DCMAKE_INSTALL_RPATH="$CONDA_PREFIX/lib" -DCMAKE_FIND_ROOT_PATH="$CONDA_PREFIX" -DLIBRARY_INSTALL_DIR="$CONDA_PREFIX/lib" $CMAKE_EXTRA_OPTS
	make -j $NUMCORES
	make install
done

echo "Building inversion modules"
cd "$D/GRALE2/inversion_modules"
if ! [ -e "build" ] ; then
	mkdir build
fi
cd build
cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" -DCMAKE_INSTALL_RPATH="$CONDA_PREFIX/lib" -DCMAKE_FIND_ROOT_PATH="$CONDA_PREFIX" -DLIBRARY_INSTALL_DIR="$CONDA_PREFIX/lib" $CMAKE_EXTRA_OPTS
make -j $NUMCORES
make install

cd "$D/GRALE2/pygrale"
CXXFLAGS="-O3 -std=c++11" ./setup.py build install

