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
import subprocess

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

if "INCLUDES" in os.environ:
    extraIncludes += os.environ["INCLUDES"].split(":")

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
        [ "grale/contourfinder.pyx", "grale/multicontourfinder.cpp" ],
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
        pyName = os.path.join(*pathParts, outName)
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
        cwd = os.getcwd()
        try:
            os.chdir(os.path.join("grale","editor","cppqt"))
            if not os.path.exists("Makefile"):
                subprocess.check_call( [ sys.executable, "configure.py" ])
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


