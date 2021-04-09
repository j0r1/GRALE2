#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import time
import subprocess

def usage(full = False):
    print(r"""Usage:
    python BuildInCondaEnvironment.py [-builddir] path
or
    python BuildInCondaEnvironment.py --help

This script attempts to compile PyGrale in an existing conda environment.

On Windows, you'll need to install the Visual Studio 2017 build tools from 
https://www.visualstudio.com/downloads/#build-tools-for-visual-studio-2017
(see also
https://stackoverflow.com/questions/57795314/are-visual-studio-2017-build-tools-still-available-for-download)
Then, in the start menu select "x64 Native Tools Command Prompt for VS2017",
and enable the conda environment by typing c:\pathtoanaconda\Scripts\activate
That way, both the build tools will be available, as well as anaconda. Once
everything has been built, you can just start the conda prompt directly.

Start with creating a new conda environment, and activate it, e.g:

    conda create -n mypygrale -c conda-forge -y python=3.8
    conda activate mypygrale

Note that we're using the conda-forge channel for packages.

Then, running the script should install/build everything in this environment,
downloading and compiling the source code in the specified path, e.g:

    python BuildInCondaEnvironment.py path_to_dir

Of course: use at your own risk!
""")

    if full:
        print("""
If '-builddir' is not specified, the path will be used as a parent directory
in which the GRALE2 repository will be cloned. This is probably what you
want.

If '-builddir' is specified, the path is interpreted as an existing empty
subdirectory of a GRALE2 repository clone, so that one directory above this
contains the main CMakeLists.txt. No git update will be performed. No other
sources will be built. This can be helpful for development purposes.

Environment variables that are used:

 - NUMCORES: if set, 'make' will be started to use this amount of
   cores for the build (with '-j' option). If not set, the number
   of cores will be detected using the python multiprocessing module

 - NOENVCHECK: don't check the conda environment

 - NOPULL: skip git pull to update used repositories

 - NOCONDAINSTALL: assume all needed dependencies are available,
   don't run 'conda install' to add any packages.

 - NOGSL: if you don't want to install a GSL package using anaconda,
   e.g. if you want to use a different one, you can set this.

 - NOMPI: if you don't want to install the openmpi package using
   anaconda, you can set this.

 - CMAKE_EXTRA_OPTS: extra options that are passed to each 'cmake'
   command. I use it to set 'ADDITIONAL_LIBRARIES_CORE' for GRALE2
   to be able to use the correct GSL libraries, with something like

    CMAKE_EXTRA_OPTS=-DADDITIONAL_LIBRARIES_CORE=-L\`gsl-config --prefix\`/lib

   (a newer GSL library needs to be used, enabled using a 'module load'
   command, but the system-wide one is used in the linking process. Seems
   similar as in e.g. https://cmake.org/pipermail/cmake/2011-June/044790.html,
   where the library's path is removed from the link command because it
   is also set in \$LIBRARY_PATH)

 - NOEMPTYDIRCHECK: in case '-builddir' is used, the specified directory
   is expected to be empty. You can set this environment variable to
   skip this check.

""")
    sys.exit(-1)

def getInstalledCondaPackages():
    return [ l.split()[0] for l in subprocess.check_output("conda list", shell=True).decode().splitlines() if not l.startswith("#") ]

def checkCall(cmd, errString):
    try:
        subprocess.check_call(cmd, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL, shell=True)
    except Exception as e:
        raise Exception(errString + " ({})".format(e))

def getUnixGeneratorAndMakeCommand():
    import multiprocessing
    return ("Unix Makefiles", [ "make", "-j", "{}".format(multiprocessing.cpu_count())], "make --version")

def getWindowsGeneratorAndMakeCommand():
    return ("NMake Makefiles", [ "nmake" ], "nmake /HELP")

getGeneratorAndMakeCommand = getWindowsGeneratorAndMakeCommand if os.name == 'nt' else getUnixGeneratorAndMakeCommand

def main():
    try:
        if len(sys.argv) != 2 and len(sys.argv) != 3:
            raise Exception("Invalid number of arguments")

        if sys.argv[1] in [ "-h", "-help", "--h", "--help" ]:
            usage(True)

        D = None
        buildDir = False
        if sys.argv[1] == "-builddir":
            if len(sys.argv) != 3:
                raise Exception("Argument -builddir specified, but no path given")
            D = sys.argv[2]
            buildDir = True
        else:
            if len(sys.argv) == 3:
                raise Exception("Too many arguments")
            D = sys.argv[1]
    except Exception as e:
        print("ERROR:", e)
        usage()

    if not D:
        raise Exception("Specified path is an empty string")
 
    generator, makeCmd, makeTestCommand = getGeneratorAndMakeCommand()

    checkCall("git --version", "Git does not appear to be available")
    checkCall("conda --version", "Need active conda environment - 'conda' command is not available")
    checkCall(makeTestCommand, "Need make/nmake to be able to build the source code")

    if not "CONDA_PREFIX" in os.environ:
        raise Exception("CONDA_PREFIX environment variable is not set, is conda environment active?")

    startDir = os.getcwd()
    if buildDir:
        if not os.path.exists(D):
            raise Exception("Specified path '{}' does not exist".format(D))

        os.chdir(D)
        if not "NOEMPTYDIRCHECK" in os.environ:
            if os.listdir():
                raise Exception("Specfified build directory is not empty")

    else:
        if not os.path.exists(D):
            os.mkdir(D)
        os.chdir(D)

    D = os.getcwd()

    if not "NOENVCHECK" in os.environ:
        if not "CONDA_DEFAULT_ENV" in os.environ:
            raise Exception("No conda default environment could be determined")

        if os.environ["CONDA_DEFAULT_ENV"] == "base":
            raise Exception("Refusing to install in base conda environment - set NOENVCHECK to override")

        print()
        print("Will install in conda environment '{}'".format(os.environ["CONDA_DEFAULT_ENV"]))
        print()
        print("Continuing in 2 seconds...")
        time.sleep(1)
        print("Continuing in 1 seconds...")
        time.sleep(1)

    startPacks = getInstalledCondaPackages()
    if not "python" in startPacks:
        print("""
Please make sure that 'python' is installed in the environment, from conda-forge. You
can run e.g.

    conda install -c conda-forge python

or

    conda install -c conda-forge python=3.8

""")
        sys.exit(-1)

    if not "NOCONDAINSTALL" in os.environ:
        condaPacks = [ "ipython", "jupyter", "astropy", "pyqt5-sip", "pyqt", "cython", "numpy", "scipy", "matplotlib",
                       "shapely", "pyopengl", "ipywidgets", "cmake" ]

        if os.name != "nt":
            condaPacks.append("compilers")

        if not "NOSGL" in os.environ:
            condaPacks.append("gsl")
        if not "NOMPI" in os.environ:
            if os.name != "nt":
                condaPacks.append("openmpi")

        remainingCondaPacks = [ ]
        for i in condaPacks:
            if not i in startPacks:
                remainingCondaPacks.append(i)

        condaPacks = remainingCondaPacks
        if not condaPacks:
            print("All required conda packages seem to be installed, continuing...")
        else:
            print("Installing {} from conda-forge".format(" ".join(condaPacks)))
            
            cmd = [ "conda", "install", "-y", "-c", "conda-forge" ] + condaPacks
            subprocess.check_call(" ".join(cmd), shell=True)

            newPacks = getInstalledCondaPackages()
            if not "compilers" in startPacks and "compilers" in newPacks:
                print("""
The conda-forge compilers package appears to be newly installed. Please exit 
and reactivate this conda environment so that the compiler will be detected
correctly in the subsequent build.

Run:

    conda deactivate
    conda activate {}

and re-run this script.

""".format(os.environ["CONDA_DEFAULT_ENV"]))
                sys.exit(-1)


    extraOpts = [] if not "CMAKE_EXTRA_OPTS" in os.environ else os.environ["CMAKE_EXTRA_OPTS"].split()

    CP = os.environ["CONDA_PREFIX"]
    if os.name == "nt":
        CP += "/Library"

    cmakeCmd = [ "cmake", "..", "-DCMAKE_BUILD_TYPE=release", "-DCMAKE_INSTALL_PREFIX=" + CP,
                 "-DCMAKE_INSTALL_RPATH=" + CP + "/lib", "-DCMAKE_FIND_ROOT_PATH=" + CP,
                 "-DLIBRARY_INSTALL_DIR=" + CP + "/lib", "-G", generator ] + extraOpts

    def configureAndInstall():
        subprocess.check_call(cmakeCmd)
        subprocess.check_call(makeCmd)
        subprocess.check_call(makeCmd + [ "install"])

    if not buildDir:

        for p in [ "ErrUt", "SerUt", "EATk", "GRALE2" ]:
            os.chdir(D)
            if not os.path.exists(p):
                print("Cloning {}".format(p))
                subprocess.check_call([ "git", "clone", "https://github.com/j0r1/" + p ])
                subprocess.check_call([ "git", "config", "pull.ff", "only" ])
                os.chdir(p)
            else:
                print("Updating {}".format(p))
                os.chdir(p)
                if not "NOPULL" in os.environ:
                    subprocess.check_call( [ "git", "pull" ])
                else:
                    print("Environment variable NOPULL was set, skipping git pull")

            print("Building " + p)
            if not os.path.exists("build"):
                os.mkdir("build")
            os.chdir("build")
            configureAndInstall()

        else:
            configureAndInstall()

    print("Building inversion modules")
    os.chdir(os.path.join("..", "inversion_modules"))
    if not os.path.exists("build"):
        os.mkdir("build")
    os.chdir("build")

    configureAndInstall()

    os.chdir(os.path.join("..", "..", "pygrale"))
    subprocess.check_call(["python", "setup.py", "build", "install"])

if __name__ == "__main__":
    main()

