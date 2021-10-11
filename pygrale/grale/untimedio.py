#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import time
from timeit import default_timer as timer
from .debuglog import debugLog

class IOException(Exception):
    pass

class IOTimeoutException(Exception):
    pass

class IOConnectionClosedException(IOException):
    def __str__(self):
        return "Connection closed"

# Set up functions B and S to convert a string to bytes or vice versa. Useful
# for both python 2 and 3 support
if sys.version_info[0] == 2:
    B = lambda s: s
    S = B
else:
    B = lambda s: bytes(s, 'UTF-8')
    S = lambda b: b.decode(encoding='UTF-8')

# A helper class to read from/write to file descriptors. Can use a timeout
class IO(object):
    def __init__(self, readFileDesc, writeFileDesc, defaultReadTimeout = 20):
        self.readFileDesc = readFileDesc
        self.writeFileDesc = writeFileDesc
        self.defaultReadTimeout = defaultReadTimeout
        self.inputClosed = False

    def writeBytes(self, b):
        if self.writeFileDesc is None:
            raise IOException("IO.writeBytes: no write file descriptor has been set")

        if os.write(self.writeFileDesc, b) != len(b):
            raise IOException("Unable to write specified number of bytes")
        
        debugLog("Writing '%s' on PID %d\n" % (repr(b), os.getpid()))

    def writeLine(self, l):
        self.writeBytes(B(l) + b"\n")

    def readBytesUntimed(self, numBytes):
        data = b""
        while len(data) < numBytes:
            data += os.read(self.readFileDesc, numBytes-len(data))
        return data

    def readBytes(self, num, timeout=-1):
        
        if self.inputClosed:
            raise IOConnectionClosedException

        if self.readFileDesc is None:
            raise IOException("IO.readBytes: no read file descriptor has been set")

        if num <= 0:
            raise IOException("IO.readBytes: invalid number of bytes")

        if timeout < 0:
            timeout = self.defaultReadTimeout

        startTime = timer()
        dt = timeout*1.0

        b = b"" 
        while len(b) < num and dt >= 0:

            x = os.read(self.readFileDesc, num-len(b))
            if len(x) == 0:
                self.inputClosed = True
                raise IOConnectionClosedException

            b += x

            dt = timeout - (timer()-startTime)

        if len(b) != num:
            raise IOTimeoutException

        debugLog("Read bytes '%s'" % repr(b))
        return b

    def readLine(self, timeout=-1):
        
        if self.inputClosed:
            raise IOConnectionClosedException

        if self.readFileDesc is None:
            raise IOException("IO.readLine: no read file descriptor has been set")

        if timeout < 0:
            timeout = self.defaultReadTimeout

        startTime = timer()
        dt = timeout*1.0

        l = b"" 
        while dt >= 0:

            c = os.read(self.readFileDesc, 1)
            if len(c) == 0:
                self.inputClosed = True
                raise IOConnectionClosedException

            if c == b"\n":
                break

            l += c

            dt = timeout - (timer()-startTime)

        l = S(l)
        debugLog("Read line '%s'" % l)
        return l

