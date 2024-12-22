#!/bin/bash -e

usage() {
cat << EOF

This is a script that attempts to create a basic PyGrale installation
in a certain directory. It downloads many components, but not GSL or
an MPI implementation.

Of course: use at your own risk

Environment variables that are used:

 - NUMCORES: if set, 'make' will be started to use this amount of
   cores for the build (with '-j' option). If not set, the number
   of cores will be detected using the python multiprocessing module

 - NOPULL: skip git pull to update used repositories

 - CMAKE_EXTRA_OPTS: extra options that are passed to each 'cmake'
   command. I use it to set 'ADDITIONAL_LIBRARIES_CORE' for GRALE2
   to be able to use the correct GSL libraries, with something like

    CMAKE_EXTRA_OPTS=-DADDITIONAL_LIBRARIES_CORE=-L\`gsl-config --prefix\`/lib

   (a newer GSL library needs to be used, enabled using a 'module load'
   command, but the system-wide one is used in the linking process. Seems
   similar as in e.g. https://cmake.org/pipermail/cmake/2011-June/044790.html,
   where the library's path is removed from the link command because it
   is also set in \$LIBRARY_PATH)

 - NOVENV: if set, no virtual environment will be used, everything will
   be installed globally in python

It may also be useful to set CC and CXX (C and C++) compiler environment
variables: cmake might detect the wrong one if more are available. I.e.
on one cluster I need to set

    export CC=icc
    export CXX=icc

to use the Intel compiler, while on another I need to set

    export CC=gcc
    export CXX=g++

to use the correct GNU compiler enabled with a 'module load' command
(cmake seems to prefer e.g. /usr/bin/cc otherwise).

This script no longer tries to build grale_editor. If you need this,
the conda build script will probably be easier to use.
EOF
}

if [ "$#" != 1 ] ; then
	echo "Please specify a virtualenv directory"
	usage
	exit -1
fi

if ! [ -z "$PYTHON" ] ; then
	echo "Specified python executable $PYTHON"
	X=`"$PYTHON" --version 2>&1 |cut -f 2 -d  " "|cut -f 1 -d .`
	if [ "$X" != 3 ] ; then
		echo "This does not seem to be a python 3 executable, exiting"
		exit -1
	fi
else
	X=`python --version 2>&1 |cut -f 2 -d  " "|cut -f 1 -d .`
	if [ "$X" = "3" ] ; then
		PYTHON=python
	else
		X=`python3 --version 2>&1 |cut -f 2 -d  " "|cut -f 1 -d .`
		if [ "$X" = "3" ] ; then
			PYTHON=python3
		else
			echo "Unable to find python 3 executable, you can try setting PYTHON environment variable"
			exit -1
		fi
	fi
	echo "Using python executable '$PYTHON'"
fi

if ! git --version ; then
	echo "Need git"
	exit -1
fi

if ! cmake --version ; then
	echo "Need cmake"
	exit -1
fi

STARTDIR=`pwd`
D="$1"
if ! [ -e "$D" ] ; then
	mkdir "$D"
fi
cd "$D"
if ! [ -e lib ] ; then
	mkdir lib
fi
if ! [ -e lib64 ] ; then
	ln -s lib lib64
fi
if ! [ -e pipcache ] ; then
	mkdir pipcache
fi

PREFIX=`pwd`
if [ -z "$NUMCORES" ] ; then
	NUMCORES=`"$PYTHON" -c "import multiprocessing;print(multiprocessing.cpu_count())"`
fi

export PATH="$PREFIX/bin:$PATH"

if [ -z "$NOVENV" ] ; then

	if ! [ -e "$PREFIX/bin/activate" ] ; then
		"$PYTHON" -m venv --copies "$PREFIX"
	fi
	source "$PREFIX/bin/activate"

	if ! [ -e "$PREFIX/bin/pip" ] ; then
		echo "Can't find pip in virtual environment"
		exit -1
	fi
	PIP=pip
else
	echo "Not activating virtual environment, detecting pip"
	if pip3 --version ; then
		PIP=pip3
	elif pip --version ; then
		PIP=pip
	else
		echo "Couldn't detect pip, aborting"
		exit -1
	fi
fi

$PIP install --cache-dir "$D/pipcache" setuptools numpy scipy astropy shapely cython matplotlib qpsolvers

if ! $PIP install --cache-dir "$D/pipcache" pycairo ; then
	echo "WARNING: pycairo could not be installed - perhaps some other package/module"
	echo "         needs to be activated first"
	echo
	echo "         You can use GRALE without this, will continue with installation in"
	echo "         a moment"
	sleep 5
fi

if ! [ -e "$PREFIX/src" ] ; then
	mkdir "$PREFIX/src"
fi

SITEPACKAGES=`echo "$PREFIX"/lib/python*/site-packages/`
echo "SITEPACKAGES=$SITEPACKAGES"

for p in ErrUt SerUt EATk GRALE2 ; do
	cd "$PREFIX/src"
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
	cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_FIND_ROOT_PATH="$PREFIX" -DCMAKE_INSTALL_RPATH="$PREFIX/lib" -DLIBRARY_INSTALL_DIR="$PREFIX/lib" $CMAKE_EXTRA_OPTS
	make -j $NUMCORES
	make install
done

cd "$PREFIX/src/GRALE2/pygrale"

if ! [ -e "Makefile" ] ; then
	./configure.py
fi

CXXFLAGS="-O3 -std=c++11" make
CXXFLAGS="-O3 -std=c++11" make install

cd "$PREFIX/src"
if ! [ -e triangle ] ; then
	mkdir triangle
fi
echo "Downloading 'triangle' sources"
cd triangle
if ! [ -e "triangle.zip" ] ; then
	curl -o triangle.zip https://netlib.org/voronoi/triangle.zip
fi

unzip -o triangle.zip

echo "Building 'triangle'"
gcc -O2 -DNO_TIMER -DLINUX -o../../bin/triangle triangle.c -lm

cat > "$PREFIX/bin/activategrale" << EOF
export PATH="$PREFIX/bin:\$PATH"
if [ -z "\$LD_LIBRARY_PATH" ] ; then
	export LD_LIBRARY_PATH="$PREFIX/lib"
else
	export LD_LIBRARY_PATH="$PREFIX/lib:\$LD_LIBRARY_PATH"
fi
EOF

if [ -z "$NOVENV" ] ; then
	cat >> "$PREFIX/bin/activategrale" << EOF
source "$PREFIX/bin/activate"
EOF
fi

cat << EOF

---------------------------------------------------------------------

Type

  source "$PREFIX/bin/activategrale"

to enable pygrale environment

EOF

