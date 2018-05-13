import sys
import os
import time

def log_dummy(*args):
    pass

def log_enabled(*args):
    t = time.time()
    f = int((t - int(t))*1000000.0)
    dateStr = time.strftime("%Y-%m-%d %H:%M:%S.{:06d} %Z: ".format(f), time.localtime(t))
    emptyStr = " "*len(dateStr)
    allText = " ".join([ "{}".format(i) for i in args ])

    prefix = dateStr
    for l in allText.split("\n"):
        sys.stderr.write(prefix + l.strip() + "\n")
        prefix = emptyStr

log = log_dummy if not "GRALESHELL_DEBUG" in os.environ else log_enabled

if __name__ == "__main__":
    log("Hello", "this", "is", "\n", 10, 12, 13)
    log("Another entry")


