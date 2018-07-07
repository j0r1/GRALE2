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
	virtualenv --always-copy "$PREFIX"
fi
source "$PREFIX/bin/activate"

$PIP install numpy scipy astropy shapely cython matplotlib

if ! [ -e "$PREFIX/src" ] ; then
	mkdir "$PREFIX/src"
fi

SITEPACKAGES=`echo "$PREFIX"/lib/python*/site-packages/`
echo "SITEPACKAGES=$SITEPACKAGES"

for p in ErrUt SerUt ENUt MOGAL GRALE2 ; do
	cd "$PREFIX/src"
	if ! [ -e $p ] ; then
		echo "Cloning $p"
		git clone https://github.com/j0r1/$p
		cd $p
	else
		echo "Updating $p"
		cd $p
		git pull
	fi

	echo "Building $p"
	if ! [ -e "build" ] ; then
		mkdir build
	fi
	cd build
	cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_FIND_ROOT_PATH="$PREFIX" -DLIBRARY_INSTALL_DIR="$PREFIX/lib"
	make -j $NUMCORES
	make install
done

echo "Building inversion modules"
cd "$PREFIX/src/GRALE2/inversion_modules"
if ! [ -e "build" ] ; then
	mkdir build
fi
cd build
cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_FIND_ROOT_PATH="$PREFIX" -DLIBRARY_INSTALL_DIR="$PREFIX/lib"
make -j $NUMCORES
make install

BUILDPYQT="no"
PYVER=`python --version 2>&1|cut -f 2 -d " "|cut -f 1 -d .`
if [ "$PYVER" = "3" ] ; then
	X=`which qmake|cat`
	if ! [ -z "$X" ] ; then
		X=`qmake --version | grep "Using Qt version 5" | cat`
		if ! [ -z "$X" ] ; then
			BUILDPYQT="yes"
			QTBASE=`which qmake|cat|rev|cut -f 3- -d /|rev`
			echo "Qt5 detected, will build PyQt5"
			echo "QTBASE=$QTBASE"
		else
			echo "Not building PyQt5, qmake found, but is not for Qt5"
		fi
	else
		echo "Not building PyQt5, qmake not found in PATH"
	fi
else
	echo "Not building PyQt5, not python 3"
fi

function download_and_extract {
	cd "$PREFIX/src"

	URL="$1"
	echo "Downloading $URL"
	FILENAME=`echo "$URL" | rev |cut -f 1 -d / | rev`
	echo "Filename is $FILENAME"
	DIRNAME="${FILENAME::${#FILENAME}-7}"
	echo "Dirname is $DIRNAME"

	if ! [ -e "$DIRNAME" ] ; then
		if ! [ -e "$FILENAME" ] ; then
			FOUND="no"
			for k in {1..10} ; do
				if curl -L -o "${FILENAME}.tmp" "$URL" ; then
					FOUND="yes"
					break
				else
					echo "Download failed, retrying in 2 seconds"
					sleep 2
				fi
			done

			if [ "$FOUND" != "yes" ] ; then
				echo "Unable to download $URL"
				exit -1
			fi

			mv "${FILENAME}.tmp" "$FILENAME"
		else
			echo "$FILENAME exists, using this" 
		fi

		tar xfz "$FILENAME"
	else
		echo "$DIRNAME exists, using this"
	fi
	
	cd "$DIRNAME"
}

if [ "$BUILDPYQT" = "yes" ] ; then
	QTVERSION=`qmake --version |cut -f 4 -d " "`
	PYQTVERSIONS=`curl https://sourceforge.net/projects/pyqt/files/PyQt5/ |grep "tr title" | grep PyQt-5 | cut -f 2 -d \" |cut -f 2 -d "-"`
	echo "Qt version detected: $QTVERSION"
	#echo "PyQt versions available for download: $PYQTVERSIONS" 

	VPYQT=$(python -c "exec(\"import sys\ncompatibleVersion = lambda x, y: '.'.join(x.split('.')[:2]) == '.'.join(y.split('.')[:2])\nqtVersion = sys.argv[1]\nfor pyqtVersion in sys.argv[2:]:\n    if compatibleVersion(qtVersion, pyqtVersion):\n        print(pyqtVersion)\n        sys.exit(0)\nraise Exception('No compatible version found')\n\")" $QTVERSION $PYQTVERSIONS )
	echo "Compatible PyQt version: $VPYQT"
	PYQTURL="https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-${VPYQT}/PyQt5_gpl-${VPYQT}.tar.gz"
	echo "PyQt URL: $PYQTURL"

	# First download PyQt so that we can check which SIP version to use
	# (I tried using the latest (4.19.10) with PyQt 5.9.2, but that caused
	# a segfault due to a sip option '-n' not being specified
	download_and_extract "$PYQTURL"
	SIPVERSION=`grep "SIP_MIN_VERSION" configure.py | head -n 1 | cut -f 2 -d "'"`
	NEEDPYQT5SIP=`grep "from PyQt5 import sip" configure.py |cat`
	echo "Needed SIP version: $SIPVERSION"

	SIPURL="https://sourceforge.net/projects/pyqt/files/sip/sip-${SIPVERSION}/sip-${SIPVERSION}.tar.gz"
	echo "SIP URL: $SIPURL"

	# For now we'll always keep this legacy system around
	# (don't yet know how the new system will need to be used)
	download_and_extract "$SIPURL"
	# avoid rebuilding if not necessary
	if ! [ -e ".buildcomplete" ] ; then
		rm -f Makefile
		python configure.py --sysroot="$PREFIX"
		make clean
		make install

		if ! [ -z "$NEEDPYQT5SIP" ] ; then
			rm -f Makefile
			python configure.py --sysroot="$PREFIX" --sip-module=PyQt5.sip
			make clean
			make install
		fi
		touch ".buildcomplete"
	fi

	download_and_extract "$PYQTURL"
	# avoid rebuilding if not necessary
	if ! [ -e ".buildcomplete" ] ; then
		if ! [ -e Makefile ] ; then
			python configure.py --confirm-license `for i in QtHelp QtMultimedia QtMultimediaWidgets QtNetwork QtPrintSupport QtQml QtQuick QtSql QtSvg QtTest QtWebKit QtWebKitWidgets QtXml QtXmlPatterns QtDesigner QAxContainer QtDBus QtSensors QtSerialPort QtX11Extras QtBluetooth QtMacExtras QtPositioning QtWinExtras QtQuickWidgets QtWebSockets QtWebChannel QtLocation QtNfc QtWebEngineCore QtWebEngine QtWebEngineWidgets Enginio; do echo "--disable $i" ; done`
			make clean
		fi
		# building with multiple cores sometimes appears to fail
		#make -j $NUMCORES install
		make install
		touch ".buildcomplete"
	fi
fi

cd "$PREFIX/src/GRALE2/pygrale"
export SIPINCLUDES="$PREFIX/share/sip/PyQt5/" 
export QT5INCLUDES="$QTBASE/include:$QTBASE/include/QtCore:$QTBASE/include/QtGui:$QTBASE/include/QtWidgets"
export QT5LIBDIRS="$QTBASE/lib"
CXXFLAGS="-O3 -std=c++11" ./setup.py build install --prefix="$PREFIX"

cd "$PREFIX/src"
if ! [ -e triangle ] ; then
	mkdir triangle
fi
echo "Downloading 'triangle' sources"
cd triangle
if ! [ -e "triangle.zip" ] ; then
	curl -o triangle.zip http://www.netlib.org/voronoi/triangle.zip
fi

unzip -o triangle.zip
if [ -z "$CC" ] ; then
	CC=gcc
fi

echo "Building 'triangle'"
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

---------------------------------------------------------------------

Type

  source "$PREFIX/bin/activategrale"

to enable pygrale environment

EOF

