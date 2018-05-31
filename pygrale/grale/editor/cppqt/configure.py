#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import subprocess
import sipconfig
from PyQt5 import QtCore
import pprint
import platform

prefix = None
for a in sys.argv[1:]:
    prefArg = "--prefix="
    if a.startswith(prefArg):
        prefix = a[len(prefArg):]
    else:
        raise Exception("Unknown argument {}".format(a))

build_file = "cppqt.sbf"

config = sipconfig.Configuration()
if prefix:
    print("Prefix set to", prefix)
    sysRoot = None
    for i in config.sip_config_args.split():
        sysRootArg="--sysroot="
        if i.startswith(sysRootArg):
            sysRoot=i[len(sysRootArg):]
    if not sysRoot:
        raise Exception("Prefix set, but sysroot could not be detected")
    
    print("Before:", config.default_mod_dir)
    config.default_mod_dir = config.default_mod_dir.replace(sysRoot, prefix)
    print("After:", config.default_mod_dir)

config.default_mod_dir = os.path.join(config.default_mod_dir, "grale", "editor", "cppqt")

macros = config.build_macros()
for x in [ "CXX", "CC" ]:
    if x in os.environ:
        macros[x] = os.environ[x]

if "CC" in os.environ:
    x = os.environ["CC"]
    macros["LINK"] = x
    macros["LINK_SHLIB"] = x

for x in [ "CXXFLAGS", "CFLAGS" ]:
    if x in os.environ:
        macros[x] += " " + os.environ[x]

config.set_build_macros(macros)
#pprint.pprint(config.build_macros())

extraCppIncludes = [ ]
extraLibraryDirs = [ ]
extraSipIncludes = [ ]
if "CONDA_PREFIX" in os.environ:
    for i in [ ".", "QtCore", "QtWidgets", "QtGui" ]:
        if platform.system() == "Windows":
            extraCppIncludes.append(os.path.join(os.environ["CONDA_PREFIX"],"Library","include","qt",i))
        else:
            extraCppIncludes.append(os.path.join(os.environ["CONDA_PREFIX"],"include","qt",i))

    if platform.system() == "Windows":
        extraLibraryDirs.append(os.path.join(os.environ["CONDA_PREFIX"], "Library", "lib"))
        extraSipIncludes.append("-I" + os.path.join(os.environ["CONDA_PREFIX"], "sip", "PyQt5"))
    else:
        extraLibraryDirs.append(os.path.join(os.environ["CONDA_PREFIX"], "lib"))
        extraSipIncludes.append("-I" + os.path.join(os.environ["CONDA_PREFIX"], "share", "sip", "PyQt5"))

if "SIPINCLUDES" in os.environ:
    extraSipIncludes += [ "-I" + p for p in os.environ["SIPINCLUDES"].split(":") ]
if "QT5LIBDIRS" in os.environ:
    extraLibraryDirs += os.environ["QT5LIBDIRS"].split(":")
if "QT5INCLUDES" in os.environ:
    extraCppIncludes += os.environ["QT5INCLUDES"].split(":")

sipFlags = QtCore.PYQT_CONFIGURATION["sip_flags"].split()
print("Generating C++ code for PyQt5 bindings")
subprocess.check_call([config.sip_bin, "-c", ".", "-b", build_file] + extraSipIncludes + sipFlags + [ "cppqt.sip"])

makefile = sipconfig.SIPModuleMakefile(config, build_file)
makefile.extra_include_dirs = extraCppIncludes
makefile.extra_lib_dirs = extraLibraryDirs
makefile.extra_libs = [ "Qt5Widgets", "Qt5Core", "Qt5Gui" ]
makefile.generate()
