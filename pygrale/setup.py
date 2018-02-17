#!/usr/bin/env python

from distutils import sysconfig
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import platform
import glob
import sys
import os
import pprint

extraFlags = [ ]
extraIncludes = [ ]
extraDefines = [ ]
libDirs = [ ]

libraries = [ "gsl", "serut", "grale2core" ]
if platform.system() == "Windows":
    libraries += [ "cblas", "errut" ]
else:
    libraries += [ "gslcblas" ]
    extraFlags += [ "-std=c++11" ]

if "CONDA_BUILD" in os.environ or "CONDA_PREFIX" in os.environ:
    prefix = None
    for n in [ "PREFIX", "CONDA_PREFIX" ]:
        if n in os.environ:
            prefix = os.environ[n]
            print("Checking prefix: " + prefix)

            l = glob.glob(os.path.join(prefix,"lib","python*","site-packages","numpy","core","include"))
            if len(l) > 0:
                extraIncludes += [ l[0] ]
            l = glob.glob(os.path.join(prefix,"lib","site-packages","numpy","core","include"))
            if len(l) > 0:
                extraIncludes += [ l[0] ]

            if extraIncludes:
                break

    if platform.system() == "Windows":
        extraIncludes += [ os.path.join(prefix, "Library", "include") ]
        libDirs += [ os.path.join(prefix, "Library", "lib") ]
    else:
        extraIncludes += [ os.path.join(prefix, "include") ]

print("Extra includes:")
pprint.pprint(extraIncludes)
print("Libraries:")
pprint.pprint(libraries)

if "CC" in os.environ and os.environ["CC"].startswith("ccache"):
    del os.environ["CC"]
if "CXX" in os.environ and os.environ["CXX"].startswith("ccache"):
    del os.environ["CXX"]

if "TEST_BUILD" in os.environ:

    pref = os.path.join(sys.exec_prefix, "pkgs/numpy-")

    def getNumPyVersion(s):                                                                             
        l = len(pref)
        m = [ x for x in map(int, s[l:].split("-")[0].split(".")) ]
        return m

    l = glob.glob(pref + "*/lib/python*/site-packages/numpy/core/include")
    l = l[sorted([ ( getNumPyVersion(l[i]), i ) for i in range(len(l)) ])[-1][-1]]

    extraIncludes += [ l, os.path.join(os.getcwd(),"../local/include/") ]

extensions = [
    Extension("grale.lenses", 
        [ "grale/lenses.pyx" ],
        include_dirs = extraIncludes,
        libraries = libraries,
        library_dirs = libDirs,
        language = "c++",
        define_macros = [ ] + extraDefines,
        extra_compile_args = extraFlags
    ),
    Extension("grale.images", 
        [ "grale/images.pyx", "grale/pylensplane.cpp" ],
        include_dirs = extraIncludes,
        libraries = libraries,
        library_dirs = libDirs,
        language = "c++",
        define_macros = [ ] + extraDefines,
        extra_compile_args = extraFlags
    ),
    Extension("grale.privutilcython", 
        [ "grale/privutilcython.pyx" ],
        include_dirs = extraIncludes,
        libraries = libraries,
        library_dirs = libDirs,
        language = "c++",
        define_macros = [ ] + extraDefines,
        extra_compile_args = extraFlags
    ),
    Extension("grale.inversionparams", 
        [ "grale/inversionparams.pyx" ],
        include_dirs = extraIncludes,
        libraries = libraries,
        library_dirs = libDirs,
        language = "c++",
        define_macros = [ ] + extraDefines,
        extra_compile_args = extraFlags
    ),
    Extension("grale.gridfunction", 
        [ "grale/gridfunction.pyx" ],
        include_dirs = extraIncludes,
        libraries = libraries,
        library_dirs = libDirs,
        language = "c++",
        define_macros = [ ] + extraDefines,
        extra_compile_args = extraFlags
    ),
    Extension("grale.contourfinder", 
        [ "grale/contourfinder.pyx" ],
        include_dirs = extraIncludes,
        libraries = libraries,
        library_dirs = libDirs,
        language = "c++",
        define_macros = [ ] + extraDefines,
        extra_compile_args = extraFlags
    ),
]

def getVersionString():
    curScript = sys.argv[0]
    scriptDir = os.path.dirname(os.path.abspath(curScript))

    dirName = os.path.basename(scriptDir)
    pyGralePrefix = "pygrale2-"
    if dirName.startswith(pyGralePrefix):
        versionString = dirName[len(pyGralePrefix):]
        # TODO: further checks?
        return versionString

    gralePrefix = "grale2-"
    parentDir = os.path.dirname(scriptDir)
    cmakeFile = os.path.join(parentDir, "CMakeLists.txt")
    if os.path.basename(parentDir).startswith(gralePrefix) and os.path.exists(cmakeFile):
        versionStr = [ l for l in open(cmakeFile).readlines() if "set(VERSION" in l ][0].replace(")","").split()[1].strip()
        return versionStr

    # Anaconda uses PKG_VERSION
    if "PKG_VERSION" in os.environ:
        return os.environ["PKG_VERSION"]

    raise Exception("Couldn't get version string for pygrale2")


versionStr = getVersionString()

pyMods = [ "grale.cosmology", "grale.plotutil", "grale.constants", "grale.renderers", "grale.timedio", "grale.untimedio",
           "grale.privutil", "grale.debuglog", "grale.feedback", "grale.bytestring", "grale.inverters",
           "grale.inversion", "grale.grid", "grale.multiplane", "grale.privimages" ]
setup(name = "grale", version = versionStr, ext_modules = cythonize(extensions), py_modules = pyMods)


