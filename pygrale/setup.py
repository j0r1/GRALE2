#!/usr/bin/env python

from __future__ import print_function
from distutils import sysconfig
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import platform
import glob
import sys
import os
import pprint
import subprocess

def addNumPyDirs(includeDirs, libDirs):
    try:
        # Try to deduce paths
        import numpy

        numpyIncludeDir = os.path.join(os.path.dirname(numpy.__file__), "core", "include")
        if os.path.exists(numpyIncludeDir):
            includeDirs.append(numpyIncludeDir)
    except Exception as e:
        print("Exception during addNumPyDirs: {}".format(e))

def addGeneralDirs(includeDirs, libDirs):
    try:
        modDir = os.path.dirname(subprocess.__file__) 
        if os.path.basename(modDir).startswith("python"):
            modDir = os.path.dirname(modDir)

        # TODO: make this work on windows as well
        prefix = os.path.dirname(modDir)
        includeDirs.append(os.path.join(prefix, "include"))
    except Exception as e:
        print("Exception during addGeneralDirs: {}".format(e))

def addGslDirs(includeDirs, libDirs):
    try:
        libs = subprocess.check_output( [ "gsl-config", "--libs"] ).decode().split()
        for l in libs:
            if l.startswith("-L"):
                libDirs.append(l[2:])

        includes = subprocess.check_output( [ "gsl-config", "--cflags" ] ).decode().split()
        for i in includes:
            if i.startswith("-I"):
                includeDirs.append(i[2:])

    except Exception as e:
        print("Exception during addGslDirs: {}".format(e))

def addGraleDirs(includeDirs, libDirs):
    try:
        graleBinDir = None
        for p in os.environ["PATH"].split(os.pathsep):
            if os.path.exists(os.path.join(p, "grale_massdens_threads")):
                graleBinDir = p
                break

        if graleBinDir:
            while graleBinDir.endswith(os.path.sep):
                graleBinDir = graleBinDir[:-1]
            graleDir = os.path.dirname(graleBinDir)

            includeDirs.append(os.path.join(graleDir, "include"))
            libDirs.append(os.path.join(graleDir, "lib"))

    except Exception as e:
        print("Exception during addGraleDirs: {}".format(e))

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

else:

    for f in [ addNumPyDirs, addGeneralDirs, addGslDirs, addGraleDirs ]:
        f(extraIncludes, libDirs)

if "INCLUDES" in os.environ:
    extraIncludes += os.environ["INCLUDES"].split(":")

print("Extra includes:")
pprint.pprint(extraIncludes)
print("Library dirs:")
pprint.pprint(libDirs)
print("Libraries:")
pprint.pprint(libraries)

if "CC" in os.environ and os.environ["CC"].startswith("ccache"):
    del os.environ["CC"]
if "CXX" in os.environ and os.environ["CXX"].startswith("ccache"):
    del os.environ["CXX"]

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
        [ "grale/contourfinder.pyx", "grale/multicontourfinder.cpp" ],
        include_dirs = extraIncludes,
        libraries = libraries,
        library_dirs = libDirs,
        language = "c++",
        define_macros = [ ] + extraDefines,
        extra_compile_args = extraFlags
    ),
    Extension("grale.multiplanecuda", 
        [ "grale/multiplanecuda.pyx" ],
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
    parentDir = os.path.dirname(scriptDir)
    cmakeFile = os.path.join(parentDir, "CMakeLists.txt")
    if os.path.exists(cmakeFile):
        versionStr = [ l for l in open(cmakeFile).readlines() if "set(VERSION" in l ][0].replace(")","").split()[1].strip()
        return versionStr

    # Anaconda uses PKG_VERSION
    if "PKG_VERSION" in os.environ:
        return os.environ["PKG_VERSION"]

    raise Exception("Couldn't get version string for pygrale2")

versionStr = getVersionString()
print("Using version string: '{}'".format(versionStr))

pyMods = [ "grale.cosmology", "grale.plotutil", "grale.constants", "grale.renderers", "grale.timedio", "grale.untimedio",
           "grale.privutil", "grale.debuglog", "grale.feedback", "grale.bytestring", "grale.inverters",
           "grale.inversion", "grale.grid", "grale.multiplane", "grale.privimages" ]
extraSetupArgs = { }

try:
    from PyQt5 import QtCore
    import sipconfig

    config = sipconfig.Configuration()
    subprocess.check_call( [ config.sip_bin, "-V" ])
    isQtAvailable = True
except Exception as e:
    print("Unable to detect PyQt5 or SIP, not building Grale editor")
    isQtAvailable = False

if isQtAvailable:
    pyMods += [
               "grale.editor.actionstack",
               "grale.editor.base",
               "grale.editor.checkqt",
               "grale.editor.contourleveldialog",
               "grale.editor.debug",
               "grale.editor.doublelineedit",
               "grale.editor.fitslayerinfodialog",
               "grale.editor.fitslistwidget",
               "grale.editor.hullprocessor",
               "grale.editor.imagelayer",
               "grale.editor.__init__",
               "grale.editor.layerlist",
               "grale.editor.listwidgetbase",
               "grale.editor.__main__",
               "grale.editor.mainwindow",
               "grale.editor.nullgriddialog",
               "grale.editor.openglhelper",
               "grale.editor.pointinfodialog",
               "grale.editor.pointslayer",
               "grale.editor.pointslistwidget",
               "grale.editor.rgblistwidget",
               "grale.editor.scenes",
               "grale.editor.tools",
               "grale.editor.exportareadialog",
               ]

    for ui in [ "grale/editor/contourleveldialog.ui",
                "grale/editor/pointinfodialog.ui",
                "grale/editor/fitslayerinfodialog.ui",
                "grale/editor/fitslistwidget.ui",
                "grale/editor/rgblistwidget.ui",
                "grale/editor/pointslistwidget.ui",
                "grale/editor/mainwindow.ui",
                "grale/editor/nullgriddialog.ui",
                "grale/editor/exportareadialog.ui",
            ]:
        pathParts = os.path.dirname(ui).split("/")
        ui = os.path.join(*ui.split("/"))
        outName = "ui_" + os.path.basename(ui)[:-3] + ".py"
        
        pyName = pathParts[0]
        for p in pathParts[1:]:
            pyName = os.path.join(pyName, p)
        pyName = os.path.join(pyName, outName)

        print(ui, "->", pyName)
        cmd = [ "-o", pyName, ui ]
        if platform.system() == "Windows":
            cmd = [ "cmd.exe", "/c", "pyuic5.bat" ] + cmd
        else:
            cmd = [ "pyuic5" ] + cmd

        if not "SKIPUIC" in os.environ:
            subprocess.check_call(cmd)

        pyMods.append("grale.editor." + os.path.basename(pyName)[:-3])

    # TODO: Windows?
    if not "scripts" in extraSetupArgs:
        extraSetupArgs["scripts"] = [ ]

    startScript = "grale/editor/grale_editor" if not platform.system() == "Windows" else "grale\\editor\\grale_editor.bat"
    extraSetupArgs["scripts"] += [ startScript ]

# Run the actual setup command
setup(name = "grale", version = versionStr, ext_modules = cythonize(extensions), py_modules = pyMods, **extraSetupArgs)

# Build and install SIP bindings for PyQt5 based Grale editor
if isQtAvailable:
    makeCmd = "nmake" if platform.system() == "Windows" else "make"
    if "build" in sys.argv or "install" in sys.argv:
        prefixArg = [ ]
        for a in sys.argv:
            if a.startswith("--prefix="):
                prefixArg.append(a)

        cwd = os.getcwd()
        try:
            os.chdir(os.path.join("grale","editor","cppqt"))
            if not os.path.exists("Makefile"):
                print("Prefix set to", prefixArg)
                subprocess.check_call( [ sys.executable, "configure.py" ] + prefixArg)
            subprocess.check_call( [ makeCmd ] )
        finally:
            os.chdir(cwd)

    if "install" in sys.argv:
        cwd = os.getcwd()
        try:
            os.chdir(os.path.join("grale","editor","cppqt"))
            subprocess.check_call( [ makeCmd, "install" ])
        finally:
            os.chdir(cwd)


