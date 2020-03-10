#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import time

# Set up a function called 'debugLog', which is enabled or disabled depending on
# the existence of the environment variable DEBUGLOG_VERBOSE
if os.environ.get("DEBUGLOG_VERBOSE") is not None:
    def debugLog(s):
        lines = s.splitlines()
        t = time.time()
        intT = int(t)
        fracT = t - intT
        pref = "%s.%03d" % (time.strftime("%H:%M:%S", time.gmtime(t)), int(fracT*1000.0))
        pref += " "
        spcs = " " * len(pref)

        sys.stderr.write(pref + lines[0] + "\n")
        for i in range(1,len(lines)):
            sys.stderr.write(spcs + lines[i] + "\n")

        sys.stderr.flush()
else:
    def debugLog(s):
        pass

