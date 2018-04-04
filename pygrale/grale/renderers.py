"""This module describes the 'renderers' you can use. By default, calculating
a lens plane mapping or a mass density on a grid, will only use a single
processor core. Functions in the :py:mod:`plotutil <grale.plotutil>` module
take a `renderer` argument, with which you can specify an external process
to be used for this purpose. Depending on the type selected, this can use
multiple threads, distributed processes using MPI, or the GPU
using OpenCL, each of which can lead to a significant speedup.

A specified renderer can either be one of the objects that are present in
this module, or it can simply be a string that describes the kind of
renderer to use. For a `lens plane renderer` these strings (case insensitive)
and objects can be:

 - "Threads" or an instance of :class:`ThreadsLensPlaneRenderer`
 - "MPI" or an instance of :py:class:`MPILensPlaneRenderer`
 - "OpenCL" or an instance of :py:class:`OpenCLLensPlaneRenderer`

For a `mass density renderer` these strings (case insensitive) and 
objects can be:

 - "Threads" or an instance of :class:`ThreadsMassDensityRenderer`
 - "MPI" or an instance of :py:class:`MPIMassDensityRenderer`

"""

from __future__ import print_function
import os
import subprocess

# For now, always use the timed version: better chance of detecting a
# crashed MPI process (since we're using a named pipe in that case, the
# crashed subprocess doesn't seem to register as a closed connection?
if not hasattr(subprocess, 'STARTUPINFO'): # Not Windows: can't use 'poll' there
    #print("Using timed IO")
    from . import timedio as timed_or_untimed_io
else:
    #print("Using untimed IO")
    from . import untimedio as timed_or_untimed_io

import struct
import time
import tempfile

debugOutput = False
debugDirectStderr = False

"""An exception of this type is raised in case an error occurs when calculating
a lens plane or mass density map using an external renderer"""
class RendererException(Exception):
    pass

class Renderer(object):
    def __init__(self, args, renderType, extraEnv = None, feedbackObject = None, rdFileDesc = None, wrFileDesc = None):
        self.renderType = renderType
        self.args = args
        #self.numX = 0
        #self.numY = 0
        #self.bottomLeft = (0,0)
        #self.topRight = (0,0)
        self.extraEnv = extraEnv
        self.version = ""
        self.feedback = feedbackObject
        self.rdFileDesc = rdFileDesc
        self.wrFileDesc = wrFileDesc

    def setFeedbackObject(self, feedbackObject):
        self.feedback = feedbackObject

    def getFeedbackObject(self):
        return self.feedback

    def render(self, lensData, bottomLeft, topRight, numX, numY):

        self.renderParams = { "bottomleft": bottomLeft, "topright": topRight, "numx": numX, "numy": numY }
        self.renderTypeIsGrid = True

        return self._renderCommon(lensData)

    def renderXYVector(self, lensData, xyArray):

        arr = xyArray.flatten().astype("double")
        if arr.shape[0] % 2 != 0:
            raise RendererException("The xyArray should contain both X and Y coordinates")

        self.renderParams = { "num": arr.shape[0]//2, "xy": arr.tobytes() }
        self.renderTypeIsGrid = False

        return self._renderCommon(lensData)

    def _renderCommon(self, lensData):

        errFile = tempfile.TemporaryFile("w+t") if not debugDirectStderr else None
        proc = None
        try:
            env = None
            if self.extraEnv is not None:
                env = { }
                for n in os.environ:
                    env[n] = os.environ[n]

                for n in self.extraEnv:
                    env[n] = self.extraEnv[n]

            #print(self.args)
            proc = subprocess.Popen(self.args, stdin = subprocess.PIPE, stdout = subprocess.PIPE, env = env, stderr = errFile)
        except Exception as e:
            raise RendererException("Unable to start renderer process (is render type supported?): {}".format(e)) 

        try:
            self.onStatus("Using rendering process for type: " + self.renderType)

            inFd = proc.stdin.fileno() if self.wrFileDesc is None else self.wrFileDesc
            outFd = proc.stdout.fileno() if self.rdFileDesc is None else self.rdFileDesc
            io = timed_or_untimed_io.IO(outFd, inFd)

            line = io.readLine(30)
            rendererId = "RENDERER:"
            if not line.startswith(rendererId):
                raise RendererException("Unexpected idenficiation from renderer process: '{}'".format(line))

            self.version = line[len(rendererId):]
            self.onStatus("Version info: " + self.version)

            io.writeLine("LENSPLANE")
            io.writeLine("LENSSIZE {}".format(len(lensData)))
            

            io.writeBytes(lensData)

            if self.renderTypeIsGrid:
                bottomLeft, topRight = self.renderParams["bottomleft"], self.renderParams["topright"]
                numX, numY = self.renderParams["numx"], self.renderParams["numy"]
                renderInfo = struct.pack("=ddddii", bottomLeft[0], bottomLeft[1], topRight[0], topRight[1], 
                                                    numX, numY)

                io.writeLine("RENDERTYPE GRID 0")
                io.writeBytes(renderInfo)
            else:
                io.writeLine("RENDERTYPE POINTVECTOR {}".format(self.renderParams["num"]))
                io.writeBytes(self.renderParams["xy"])
            
            statusStr = "STATUS:"
            progressStr = "PROGRESS:"
            resultStr = "RESULT:"
            while True:
                # Is 600 secs long enough?
                # TODO: add a keepalive message to be able to reduce this
                line = io.readLine(600) 
                if line.startswith(statusStr):
                    s = line[len(statusStr):]
                    self.onStatus(s)
                elif line.startswith(progressStr):
                    [ cur, target ] = map(float, line[len(progressStr):].split())
                    x = min(max(int(round(100.0*cur/target)),0),100)
                    self.onProgress(x)
                elif line.startswith(resultStr):
                    numBytes = int(line[len(resultStr):])
                    break
                else:
                    raise RendererException("Unexpected line '{}'".format(line))
                
            # Read the actual rendered data
            renderData = b""
            while len(renderData) < numBytes:
                renderData += os.read(outFd, numBytes-len(renderData))

            #print("len(renderData) =", len(renderData))
            
            io.writeLine("EXIT")
            #time.sleep(10)

            #print("Waiting...")
            proc.wait()
            #print("Done waiting")

            #self.numX = numX
            #self.numY = numY
            #self.bottomLeft = self.bottomLeft
            #self.topRight = self.topRight
        except RendererException:
            raise
        except Exception as e:
            errFile.flush()
            errFile.seek(0)
            errData = errFile.read()
            errLines = [ l.strip() for l in errData.splitlines() if len(l.strip()) > 0 ]
            if len(errLines) > 0:
                errLine = "Last line from error log: " + errLines[-1]
            else:
                errLine = "Something went wrong, but don't know what"
            raise RendererException(errLine)
        finally:
            try:
                proc.stdin.close()
                proc.stdout.close()
            except:
                pass
            time.sleep(0.1)
            try:
                from . import privutil
                privutil.terminateProcess(proc, feedbackObject = self.feedback)
            except Exception as e:
                print("Ignoring exception when terminating program: " + str(e))
        
            if debugOutput: # for debugging
                errFile.flush()
                errFile.seek(0)
                errData = errFile.read()
                print()
                print("LOG:")
                print(errData)

        return renderData

    def onStatus(self, s):
        try:
            if self.feedback:
                self.feedback.onStatus(s)
        except Exception as e:
            print("Warning: ignoring exception ({}) in onStatus".format(e))

    def onProgress(self, x):
        try:
            if self.feedback:
                self.feedback.onProgress(x)
        except Exception as e:
            print("Warning: ignoring exception ({}) in onProgress".format(e))


class ThreadsLensPlaneRenderer(Renderer):
    """Renderer that uses multiple threads to speed up the calculation of a
    lens plane."""

    def __init__(self, numProcesses = None, feedbackObject = None):
        """If `numProcesses` is set, a specific number of threads will be requested.
        A specific :py:mod:`feedback <grale.feedback>` object can be specified for 
        status and progress updates."""

        from .inverters import _getNumHelpers

        numHelpers = _getNumHelpers(numProcesses)
        if numHelpers < 1 or numHelpers > 256:
            raise RendererException("The number of processes should be at least one, and at most 256")

        e = { "GRALE_NUMTHREADS": str(numHelpers) }
        super(ThreadsLensPlaneRenderer, self).__init__([ "grale_lensplane_threads" ], "LENSPLANE", extraEnv = e, feedbackObject = feedbackObject)

class MPILensPlaneRenderer(Renderer):
    """Renderer that uses multiple processes (using MPI) to speed up the 
    calculation of a lens plane."""

    def __init__(self, numProcesses = None, feedbackObject = None):
        """If `numProcesses` is set, a specific number of threads will be requested by
        passing the ``-np`` option to ``mpirun``. A specific :py:mod:`feedback <grale.feedback>`
        object can be specified for status and progress updates."""

        from . import privutil
        self.pipePair = privutil.PipePair() # Keep it for the lifetime of the object
        pp = self.pipePair
    
        npArgs = [ ] if numProcesses is None else [ "-np", str(numProcesses) ]
        super(MPILensPlaneRenderer, self).__init__([ "mpirun"] + npArgs + [ "grale_lensplane_mpi", pp.wrFileName, pp.rdFileName ], 
                                                   "LENSPLANE", feedbackObject = feedbackObject,
                                                   rdFileDesc = pp.rdPipeDesc, wrFileDesc = pp.wrPipeDesc)

class OpenCLLensPlaneRenderer(Renderer):
    """Renderer that uses OpenCL to use your GPU to speed up the 
    calculation of a lens plane."""
    
    def __init__(self, feedbackObject = None):
        """A specific :py:mod:`feedback <grale.feedback>` object can be specified 
        for status and progress updates."""
        
        super(OpenCLLensPlaneRenderer, self).__init__([ "grale_lensplane_opencl" ], "LENSPLANE", feedbackObject = feedbackObject)

class ThreadsMassDensityRenderer(Renderer):
    """Renderer that uses multiple threads to speed up the calculation of a
    mass density map."""

    def __init__(self, numProcesses = None, feedbackObject = None):
        """If `numProcesses` is set, a specific number of threads will be requested.
        A specific :py:mod:`feedback <grale.feedback>` object can be specified for 
        status and progress updates."""

        from .inverters import _getNumHelpers

        numHelpers = _getNumHelpers(numProcesses)
        if numHelpers < 1 or numHelpers > 256:
            raise RendererException("The number of processes should be at least one, and at most 256")

        e = { "GRALE_NUMTHREADS": str(numHelpers) }
        super(ThreadsMassDensityRenderer, self).__init__([ "grale_massdens_threads" ], "MASSDENS", extraEnv = e, feedbackObject = feedbackObject)

class MPIMassDensityRenderer(Renderer):
    """Renderer that uses multiple processes (using MPI) to speed up the 
    calculation of a mass density map."""
    
    def __init__(self, numProcesses = None, feedbackObject = None):
        """If `numProcesses` is set, a specific number of threads will be requested by
        passing the ``-np`` option to ``mpirun``. A specific :py:mod:`feedback <grale.feedback>`
        object can be specified for status and progress updates."""
        
        from . import privutil
        self.pipePair = privutil.PipePair() # Keep it for the lifetime of the object
        pp = self.pipePair

        npArgs = [ ] if numProcesses is None else [ "-np", str(numProcesses) ]
        super(MPIMassDensityRenderer, self).__init__([ "mpirun"] + npArgs + [ "grale_massdens_mpi", pp.wrFileName, pp.rdFileName ], 
                                                     "MASSDENS", feedbackObject = feedbackObject,
                                                     rdFileDesc = pp.rdPipeDesc, wrFileDesc = pp.wrPipeDesc)

_defaultMassRenderer = [ None ]
_defaultLensRenderer = [ None ]

def getDefaultLensPlaneRenderer():
    """TODO:"""
    return _defaultLensRenderer[0]

def getDefaultMassRenderer():
    """TODO:"""
    return _defaultMassRenderer[0]

def setDefaultLensPlaneRenderer(x):
    """TODO:"""
    _defaultLensRenderer[0] = x

def setDefaultMassRenderer(x):
    """TODO:"""
    _defaultMassRenderer[0] = x

