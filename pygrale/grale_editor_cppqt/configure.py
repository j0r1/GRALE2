#!/usr/bin/env python
import generate_toml
import distutils.sysconfig
import os
import site
import sys
import pprint
import platform

def generateMakefile():
    # Defaults for posix
    makeOpts = {
        "makeCmd": "make",
        "copyCmd" : "cp",
        "sep": os.sep,
        "extSuffix": distutils.sysconfig.get_config_var("EXT_SUFFIX")
    }    

    if platform.system() == "Windows":
        makeOpts["copyCmd"] = "copy"
        makeOpts["makeCmd"] = "nmake"

        if "CONDA_PREFIX" in os.environ:
            # Seems that just .pyd is used then - is there a better way to detect this
            # using the SIP tools? Is this always the case in Windows?
            makeOpts["extSuffix"] = ".pyd"

    open("Makefile", "wt").write("""

EXTSUFFIX={extSuffix}
COPY={copyCmd}
MAKE={makeCmd}

MODFILE=grale_editor_cppqt$(EXTSUFFIX)

all: $(MODFILE)

install:
	cd build && $(MAKE) install

$(MODFILE): build{sep}grale_editor_cppqt{sep}$(MODFILE)
	$(COPY) "build{sep}grale_editor_cppqt{sep}$(MODFILE)" "$(MODFILE)"

build{sep}grale_editor_cppqt{sep}$(MODFILE): pyproject.toml
	sip-build

""".format(**makeOpts))

def main():

    if not os.path.exists("pyproject.toml"):
        generate_toml.genToml(sys.argv[1])
    else:
        print("pyproject.toml already exists, not recreating")

    if not os.path.exists("Makefile"):
        generateMakefile()
    else:
        print("Makefile already exists, not recreating")

if __name__ == "__main__":
    main()
