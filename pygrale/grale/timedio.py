#!/usr/bin/env python

from __future__ import print_function
import os
import select
import sys
import time
from timeit import default_timer as timer
from .debuglog import debugLog
from .bytestring import B, S

class IOException(Exception):
    pass

class IOTimeoutException(Exception):
    pass

class IOConnectionClosedException(IOException):
    def __str__(self):
        return "Connection closed"

# A helper class to read from/write to file descriptors. Can use a timeout
class IO(object):
    def __init__(self, readFileDesc, writeFileDesc, defaultReadTimeout = 20):
        self.readFileDesc = readFileDesc
        self.writeFileDesc = writeFileDesc
        self.defaultReadTimeout = defaultReadTimeout
        self.inputClosed = False

        if readFileDesc is not None:
            self.pollObject = select.poll()
            self.pollObject.register(readFileDesc, select.POLLIN|select.POLLERR|select.POLLHUP)

    def writeBytes(self, b):
        if self.writeFileDesc is None:
            raise IOException("IO.writeBytes: no write file descriptor has been set")

        if os.write(self.writeFileDesc, b) != len(b):
            raise IOException("Unable to write specified number of bytes")
        
        debugLog("Writing '%s' on PID %d\n" % (repr(b), os.getpid()))

    def writeLine(self, l):
        self.writeBytes(B(l) + b"\n")

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

            r = self.pollObject.poll(dt*1000.0)
            if len(r) == 0:
                raise IOTimeoutException
            if r[0][1]&select.POLLIN == 0:
                self.inputClosed = True
                raise IOConnectionClosedException
        
            x = os.read(self.readFileDesc, 1)
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

        errFlags = select.POLLERR | select.POLLHUP

        l = b"" 
        while dt >= 0:

            r = self.pollObject.poll(dt*1000.0)
            if len(r) == 0:
                raise IOTimeoutException
            if r[0][1]&select.POLLIN == 0 or r[0][1]&errFlags != 0:

                # Connection is closed
                if len(l) != 0:
                    break

                self.inputClosed = True
                raise IOConnectionClosedException

            c = os.read(self.readFileDesc, 1)
            if c == b"\n":
                break

            l += c

            dt = timeout - (timer()-startTime)

        l = S(l)
        debugLog("Read line '%s'" % l)
        return l

