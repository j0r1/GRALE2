#!/usr/bin/env python

import sys

if sys.version_info.major < 3:
    print("Major python version should be at least 3, version 2 is no longer supported")
    sys.exit(-1)

import platform
import glob
import os
import subprocess
import json

def addNumPyDirs(includeDirs, libDirs):
    try:
        # Try to deduce paths
        import numpy

        for core in [ "core", "_core" ]:
            numpyIncludeDir = os.path.join(os.path.dirname(numpy.__file__), core, "include")
            if os.path.exists(numpyIncludeDir):
                includeDirs.append(numpyIncludeDir)
                break

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

def checkQtAvailable():
    try:
        from PyQt5 import QtCore
        import sipbuild

        subprocess.check_call( [ "sip-build", "-V" ])
        return True
    except Exception as e:
        print("Unable to detect PyQt5 or SIP, not building Grale editor")
        return False

def main():

    extraFlags = [ ]
    extraIncludes = [ ]
    extraDefines = [ ]
    libDirs = [ ]
    libraries = [ "gsl", "serut", "grale2core" ]

    if platform.system() == "Windows":
        #libraries += [ "gslcblas", "errut" ]
        libraries += [ "gslcblas" ]
    else:
        libraries += [ "gslcblas" ]
        extraFlags += [ "-std=c++17" ]

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
                l = glob.glob(os.path.join(prefix,"lib","python*","site-packages","numpy","_core","include"))
                if len(l) > 0:
                    extraIncludes += [ l[0] ]
                l = glob.glob(os.path.join(prefix,"lib","site-packages","numpy","_core","include"))
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

    pyMods = [ "grale.cosmology", "grale.plotutil", "grale.constants", "grale.renderers", "grale.timedio", "grale.untimedio",
               "grale.privutil", "grale.debuglog", "grale.feedback", "grale.bytestring", "grale.inverters",
               "grale.inversion", "grale.grid", "grale.multiplane", "grale.privimages", "grale.privlenses",
               "grale.util", "grale.all", "grale.all_nb", "grale.lensinfocache", "grale.paramdesc", "grale.ltparamdesc" ]

    extraSetupArgs = { "scripts": [ os.path.join("scripts", "grale_socket_to_mpi.py")]}

    versionStr = getVersionString()

    isQtAvailable = checkQtAvailable()

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
                   "grale.editor.backproject",
                   "grale.editor.backgroundprocessdialog",
                   "grale.editor.backprojretracedialog",
                   "grale.editor.helpdialog",
                   "grale.editor.backprojectwidget",
                   "grale.editor.backprojectsettingsdialog",
                   "grale.editor.imgregionsettingsdialog",
                   "grale.editor.noninteractive",
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
                    "grale/editor/backgroundprocessdialog.ui",
                    "grale/editor/backprojretracedialog.ui",
                    "grale/editor/helpdialog.ui",
                    "grale/editor/backprojectwidget.ui",
                    "grale/editor/backprojectsettingsdialog.ui",
                    "grale/editor/imgregionsettingsdialog.ui",
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

        # Build the Qt native code, create python version with SIP
        curDir = os.getcwd()
        try:
            os.chdir(os.path.join("grale_editor_cppqt"))
            makeCmd = "nmake" if platform.system() == "Windows" else "make"
            subprocess.check_call( [ sys.executable, "configure.py", versionStr])
        finally:
            os.chdir(curDir)

    # extraIncludes, libraries, libDirs, extraDefines, extraFlags, pyMods, isQtAvailable

    json.dump({
        "extraIncludes": extraIncludes,
        "libraries": libraries,
        "libDirs": libDirs,
        "extraDefines": extraDefines,
        "extraFlags": extraFlags,
        "extraSetupArgs": extraSetupArgs,
        "pyMods": pyMods,
        "versionStr": versionStr
        },open("config.json", "wt"),indent=2)

    qtBuildCmd, qtInstallCmd = "", ""
    if isQtAvailable:
        qtBuildCmd = "\tcd grale_editor_cppqt && {}".format(makeCmd)
        qtInstallCmd = "\tcd grale_editor_cppqt && {} install".format(makeCmd)

    open("Makefile", "wt").write("""

all: grale/lenses.pyi grale/images.pyi grale/gridfunction.pyi
	python setup.py build
{}

install:
	python -m pip -v install .
{}
                                 
grale/lenses.pyi: grale/lenses.pyx grale/privlenses.py
	python generate_pyi.py grale/lenses.pyx grale/privlenses.py grale/lenses.pyi
                                 
grale/images.pyi: grale/images.pyx grale/privimages.py
	python generate_pyi.py grale/images.pyx grale/privimages.py grale/images.pyi
                                 
grale/gridfunction.pyi: grale/gridfunction.pyx
	python generate_pyi.py grale/gridfunction.pyx grale/gridfunction.pyi

""".format(qtBuildCmd, qtInstallCmd))
    

if __name__ == "__main__":
    main()
