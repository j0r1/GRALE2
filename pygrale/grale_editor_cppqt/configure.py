#!/usr/bin/env python
import generate_toml
import distutils.sysconfig
import os
import site
import sys

def generateMakefile():
    extSuffix = distutils.sysconfig.get_config_var("EXT_SUFFIX")
    sitePack = site.getsitepackages()

    if not sitePack:
        raise Exception("Can't detect site-packages dir")

    if len(sitePack) > 1:
        print("WARNING: more than one site-packages dir was detected")

    # TODO: cp/copy POSIX/Windows
    # TODO: windows/linux dirsep
    # TODO: make/nmake command
    open("Makefile", "wt").write("""

EXTSUFFIX={}

MODFILE=grale_editor_cppqt$(EXTSUFFIX)

all: $(MODFILE)

install:
	make -C build install

$(MODFILE): build/grale_editor_cppqt/$(MODFILE)
	cp "build/grale_editor_cppqt/$(MODFILE)" "$(MODFILE)"

build/grale_editor_cppqt/$(MODFILE): pyproject.toml
	sip-build

""".format(extSuffix))


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
