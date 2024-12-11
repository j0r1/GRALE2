#!/usr/bin/env python

import sys

if sys.version_info.major < 3:
    print("Major python version should be at least 3, version 2 is no longer supported")
    sys.exit(-1)

import sysconfig
from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import pprint
import json
import os

if not os.path.exists("config.json"):
    print("config.json doesn't exist, run configure.py first")
    sys.exit(-1)

cfg = json.load(open("config.json", "rt"))
print("Loaded settings:")
pprint.pprint(cfg)

if "CC" in os.environ and os.environ["CC"].startswith("ccache"):
    del os.environ["CC"]
if "CXX" in os.environ and os.environ["CXX"].startswith("ccache"):
    del os.environ["CXX"]

extensions = [
    Extension("grale.lenses", 
        [ "grale/lenses.pyx", "grale/threadslenscalc.cpp" ],
        include_dirs = cfg["extraIncludes"],
        libraries = cfg["libraries"],
        library_dirs = cfg["libDirs"],
        language = "c++",
        define_macros = [ ] + cfg["extraDefines"],
        extra_compile_args = cfg["extraFlags"]
    ),
    Extension("grale.images", 
        [ "grale/images.pyx", "grale/pylensplane.cpp" ],
        include_dirs = cfg["extraIncludes"],
        libraries = cfg["libraries"],
        library_dirs = cfg["libDirs"],
        language = "c++",
        define_macros = [ ] + cfg["extraDefines"],
        extra_compile_args = cfg["extraFlags"]
    ),
    Extension("grale.inversionparams", 
        [ "grale/inversionparams.pyx" ],
        include_dirs = cfg["extraIncludes"],
        libraries = cfg["libraries"],
        library_dirs = cfg["libDirs"],
        language = "c++",
        define_macros = [ ] + cfg["extraDefines"],
        extra_compile_args = cfg["extraFlags"]
    ),
    Extension("grale.gridfunction", 
        [ "grale/gridfunction.pyx" ],
        include_dirs = cfg["extraIncludes"],
        libraries = cfg["libraries"],
        library_dirs = cfg["libDirs"],
        language = "c++",
        define_macros = [ ] + cfg["extraDefines"],
        extra_compile_args = cfg["extraFlags"]
    ),
    Extension("grale.contourfinder", 
        [ "grale/contourfinder.pyx", "grale/multicontourfinder.cpp" ],
        include_dirs = cfg["extraIncludes"],
        libraries = cfg["libraries"],
        library_dirs = cfg["libDirs"],
        language = "c++",
        define_macros = [ ] + cfg["extraDefines"],
        extra_compile_args = cfg["extraFlags"]
    ),
    Extension("grale.quadprogmatrix",
        [ "grale/quadprogmatrix.pyx", "grale/qpmatrix.cpp" ],
        include_dirs = cfg["extraIncludes"],
        libraries = cfg["libraries"],
        library_dirs = cfg["libDirs"],
        language = "c++",
        define_macros = [ ] + cfg["extraDefines"],
        extra_compile_args = cfg["extraFlags"]
    ),
]

pyMods = cfg["pyMods"]
extraSetupArgs = cfg["extraSetupArgs"]
versionStr = cfg["versionStr"]
print("Using version string: '{}'".format(versionStr))

# Run the actual setup command
setupInf = setup(name = "grale",
        version = versionStr, ext_modules = cythonize(extensions, language_level="3"), 
        py_modules = pyMods,
        packages = ["grale"],
        package_dir = {"": "." },
        package_data = {"grale": ["lenses.pyi", "images.pyi", "gridfunction.pyi"]},
        **extraSetupArgs)

