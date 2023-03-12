#!/usr/bin/env python
import sys
import site
import pprint
import os
import subprocess
import glob
import platform

template = '''
# Specify sip v6 as the build system for the package.
[build-system]
requires = ["sip >=6, <7", "PyQt-builder >=1.9, <2" ]
build-backend = "sipbuild.api"

# Specify the PEP 566 metadata for the project.
[tool.sip.metadata]
name = "grale_editor_cppqt"
version = "{}"

[tool.sip]
project-factory = "pyqtbuild:PyQtProject"

# Configure the building of the fib bindings.
[tool.sip.bindings.grale_editor_cppqt]
headers = [ "cppqt.h" ]
include-dirs = {} 
library-dirs = [ "{}" ]
libraries = [ "Qt5Widgets{}", "Qt5Core{}", "Qt5Gui{}" ]
'''

def getQtCoreTomlPath():
    qtCoreToml = []
    sitePack = site.getsitepackages()
    for p in sitePack:
        corePath = os.path.join(p, "PyQt5", "bindings", "QtCore", "QtCore.toml")
        if os.path.exists(corePath):
            qtCoreToml.append(corePath)

    if not qtCoreToml:
        raise Exception("No QtCore.toml could be located in {}".format(sitePack))
    if len(qtCoreToml) > 1:
        raise Exception("More than one QtCore.toml found in {}".format(sitePack))

    return qtCoreToml[0]

def getLibDir():
    o = subprocess.check_output(["qmake", "-v"])
    o = o.splitlines()
    # QMake version 
    assert o[0].startswith(b"QMake version "), "Unexpected output of qmake, can't detect version"
    qmakeVersion = o[0].strip().split()[2]
    print("Detected qmake version", qmakeVersion)

    # Using Qt version 
    qtVStart = b"Using Qt version "
    assert o[1].startswith(qtVStart), "Unexpected output of qmake, can't get qt version"

    vPart = o[1][len(qtVStart):]
    idx = vPart.find(b" ")
    assert idx > 0, "Can't find end of version string in '{}'".format(vPart)

    qtVersion = vPart[:idx]
    print("Qt version detected is", qtVersion)
    libPart = vPart[idx:]
    assert libPart.startswith(b" in "),"Unexpected output of qmake, can't get lib dir"

    libDir = libPart[4:].replace(b"\r", b"").replace(b"\n", b"")
    print("Lib in '{}'".format(libDir))

    extraLibSuffix = ""
    if "CONDA_PREFIX" in os.environ and platform.system() == "Windows":
        # It appears that on windows, with conda-forge, eg qt5widgets_conda.lib is used
        widgets = glob.glob(os.path.join(libDir, b"Qt5Widgets*.lib"))
        if len(widgets) == 1:
            extraLibSuffix = os.path.basename(widgets[0])[10:-4].decode()
            print("Detected extra library suffix", extraLibSuffix)
        elif len(widgets) == 0:
            print("WARNING: No Qt5Widgets lib found, assuming no extra suffix")
        else:
            print("WARNING: More than one Qt5Widgets lib found, can't deduce extra suffix")
            pprint.pprint(widgets)
    
    print(widgets)

    return qmakeVersion, qtVersion, libDir, extraLibSuffix

def getIncDirs(libDir):
    incDir = os.path.join(os.path.dirname(libDir), b"include", b"qt")
    assert os.path.exists(os.path.join(incDir, b"QtCore", b"QtCore")), "Can't locate Qt include dir based on lib dir"

    return [
        incDir,
        os.path.join(incDir, b"QtCore"),
        os.path.join(incDir, b"QtWidgets"),
        os.path.join(incDir, b"QtGui"),
    ]

def genToml(version):
    qmakeVersion, qtVersion, libDir, extraLibSuffix = getLibDir()
    incDirs = getIncDirs(libDir)
    incDirs = [ i.decode() for i in incDirs ]
    incDirs.append(os.getcwd())

    toml = template.format(version,incDirs, libDir.decode(), extraLibSuffix, extraLibSuffix, extraLibSuffix)
    open("pyproject.toml", "wt").write(toml)
    print("Wrote to pyproject.toml")

def main():
    genToml(sys.argv[1])

if __name__ == "__main__":
    main()
