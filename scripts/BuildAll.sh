#!/bin/bash -e

usage() {
cat << EOF

This is a script that attempts to create a basic PyGrale installation
in a certain directory. It downloads many components, but not GSL or
an MPI implementation.

Of course: use at your own risk

EOF
}

if [ "$#" != 1 ] ; then
	echo "Please specify a virtualenv directory"
	usage
	exit -1
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
PREFIX=`pwd`
NUMCORES=`python -c "import multiprocessing;print(multiprocessing.cpu_count())"`

export PATH="$PREFIX/bin:$PATH"

if pip3 --version ; then
	PIP=pip3
elif pip --version ; then
	PIP=pip
else
	echo "Couldn't detect pip, aborting"
	exit -1
fi

if ! [ -e "$PREFIX/bin/activate" ] ; then
	if ! virtualenv --version ; then
		echo "Virtualenv not found, installing"
		$PIP install virtualenv --prefix="$PREFIX" -I
		export PYTHONPATH=`echo $PREFIX/lib/*/site-packages/`
		echo $PYTHONPATH
	fi
	virtualenv "$PREFIX"
fi
source "$PREFIX/bin/activate"

$PIP install numpy scipy astropy shapely cython

if ! [ -e "$PREFIX/src" ] ; then
	mkdir "$PREFIX/src"
fi

for p in ErrUt SerUt ENUt MOGAL GRALE2 ; do
	cd "$PREFIX/src"
	if ! [ -e $p ] ; then
		git clone https://github.com/j0r1/$p
		cd $p
	else
		cd $p
		git pull
	fi

	if ! [ -e "build" ] ; then
		mkdir build
	fi
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_FIND_ROOT_PATH="$PREFIX" -DLIBRARY_INSTALL_DIR="$PREFIX/lib"
	make -j $NUMCORES
	make install
done

cd "$PREFIX/src/GRALE2/inversion_modules"
if ! [ -e "build" ] ; then
	mkdir build
fi
cd build
cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_FIND_ROOT_PATH="$PREFIX" -DLIBRARY_INSTALL_DIR="$PREFIX/lib"
make -j $NUMCORES
make install

cd "$PREFIX/src/GRALE2/pygrale"
./setup.py build install

cd "$PREFIX/src"
if ! [ -e triangle ] ; then
	mkdir triangle
fi
cd triangle
if ! [ -e "triangle.zip" ] ; then
	curl -o triangle.zip http://www.netlib.org/voronoi/triangle.zip
fi

unzip -o triangle.zip
if [ -z "$CC" ] ; then
	CC=gcc
fi

$CC -O2 -DNO_TIMER -DLINUX -o../../bin/triangle triangle.c -lm

cat > "$PREFIX/bin/activategrale" << EOF
source "$PREFIX/bin/activate"
export PATH="$PREFIX/bin:$PATH"
if [ -z "$LD_LIBRARY_PATH" ] ; then
	export LD_LIBRARY_PATH="$PREFIX/lib"
else
	export LD_LIBRARY_PATH="$PREFIX/lib:$LD_LIBRARY_PATH"
fi
EOF

cat << EOF
Type

  source "$PREFIX/bin/activategrale"

to enable pygrale environment

EOF

