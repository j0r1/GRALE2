#!/bin/bash -e

usage() {
cat << EOF

./BuildInCondaEnvironment.sh [-builddir] path

This script attempts to compile PyGrale in an existing conda environment.

If '-builddir' is not specified, the path will be used as a parent directory
in which the GRALE2 repository will be cloned.

If '-builddir' is specified, the path is interpreted as an existing empty
subdirectory of a GRALE2 repository clone, so that one directory above this
contains the main CMakeLists.txt. No git update will be performed.

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

 - NOEMPTYDIRCHECK: in case '-builddir' is used, the specified directory
   is expected to be empty. You can set this environment variable to
   skip this check.

EOF
}

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ] ; then
	echo "Invalid number of arguments"
	usage
	exit -1
fi

if [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--h" ] || [ "$1" == "--help" ] ; then
	usage
	exit -1
fi

if [ "$1" == "-builddir" ] ; then
	if [ "$#" != 2 ] ; then
		echo "'-builddir' specified but no path given"
		usage
		exit -1
	fi

	D="$2"
	BUILDDIR="yes"
else
	if [ "$#" = 2 ] ; then
		echo "Too many arguments"
		usage
		exit -1
	fi

	D="$1"
	BUILDDIR="no"
fi

if [ -z "$D" ] ; then
	echo "Specified path is empty"
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

if [ "$BUILDDIR" = "yes" ] ; then
	if ! [ -e "$D" ] ; then
		echo "Specified path $D does not exist"
		exit -1
	fi

	cd "$D"

	if [ -z "$NOEMPTYDIRCHECK" ] ; then
		if ! [ -z "`ls -A`" ] ; then # Check directory empty
			echo "Specified build directory is not empty"
			exit -1
		fi
	fi
else
	if ! [ -e "$D" ] ; then
		mkdir "$D"
	fi

	cd "$D"
fi
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
	if [ -z "$NOGSL" ] ; then
		OWNPACKS="$OWNPACKS gsl"
	fi
	if [ -z "$NOMPI" ] ; then
		OWNPACKS="$OWNPACKS openmpi"
	fi

	echo "OWNPACKS=$OWNPACKS"
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

if [ "$BUILDDIR" != "yes" ] ; then

	p="GRALE2"
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
fi

echo "CXX=$CXX"
# Reactivate current environment, otherwise $CC/$CXX may not get set correctly
source ${CONDA_EXE::-5}activate "$CONDA_DEFAULT_ENV"
echo "CXX=$CXX"

cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" -DCMAKE_INSTALL_RPATH="$CONDA_PREFIX/lib" -DCMAKE_FIND_ROOT_PATH="$CONDA_PREFIX" -DLIBRARY_INSTALL_DIR="$CONDA_PREFIX/lib" $CMAKE_EXTRA_OPTS
make -j $NUMCORES
make install

echo "Building inversion modules"
cd "../inversion_modules"
if ! [ -e "build" ] ; then
	mkdir build
fi
cd build
cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" -DCMAKE_INSTALL_RPATH="$CONDA_PREFIX/lib" -DCMAKE_FIND_ROOT_PATH="$CONDA_PREFIX" -DLIBRARY_INSTALL_DIR="$CONDA_PREFIX/lib" $CMAKE_EXTRA_OPTS
make -j $NUMCORES
make install

cd "../../pygrale/"
CXXFLAGS="-O3 -std=c++11" ./setup.py build install
