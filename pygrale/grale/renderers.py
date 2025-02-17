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
 - "NC:portnumber" or an instance of :py:class:`NetcatLensPlaneRenderer`

For a `mass density renderer` these strings (case insensitive) and 
objects can be:

 - "Threads" or an instance of :class:`ThreadsMassDensityRenderer`
 - "MPI" or an instance of :py:class:`MPIMassDensityRenderer`
 - "NC:portnumber" or an instance of :py:class:`NetcatMassDensityRenderer`

"""

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
from . import privutil
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
        try:
            env = privutil._mergeExtraEnvironmentVariables(self.extraEnv)
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
            errStr = "ERROR:"
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
                elif line.startswith(errStr):
                    raise RendererException(line[len(errStr):].strip())
                else:
                    raise RendererException("Unexpected line '{}'".format(line))
                
            # Read the actual rendered data
            renderData = io.readBytesUntimed(numBytes)
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
            raise RendererException(privutil._getErrorLineFromErrFile(errFile, debugOutput))
        finally:
            privutil._closeStdInOutAndterminateProcess(proc, self.feedback, debugOutput)
        
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

class NetcatRendererBase(Renderer):
    def __init__(self, port, feedbackObject, typeString):
        super(NetcatRendererBase, self).__init__([ "nc", "localhost", str(port) ], typeString, feedbackObject = feedbackObject)

class NetcatMassDensityRenderer(NetcatRendererBase):
    """Forward rendering commands to another process over a TCP connection
    using the `netcat (nc) <https://linux.die.net/man/1/nc>`_ tool.
    See :py:class:`NetcatLensPlaneRenderer` for examples."""

    def __init__(self, port, feedbackObject = None):
        """Connect to localhost (127.0.0.1) on the specified `port` number. The
        idea is to use SSH tunnels to connect to another machine."""
        super(NetcatMassDensityRenderer, self).__init__(port, feedbackObject, "MASSDENS")

class NetcatLensPlaneRenderer(NetcatRendererBase):
    """Forward rendering commands to another process over a TCP connection
    using the `netcat (nc) <https://linux.die.net/man/1/nc>`_ tool.

    On the other machine, you can use `socat <https://linux.die.net/man/1/socat>`_
    to redirect stdout/stdin of one of the render programs over the connection.
    For example

    .. code-block:: bash
    
        export GRALE_NUMTHREADS=36
        socat TCP4-LISTEN:9999,reuseaddr,fork EXEC:grale_lensplane_threads

    would wait for incoming TCP connections on port 9999, and start the thread
    based lens plane renderer when an incoming connection was detected.

    To run one of the MPI based renderers, `socat` cannot be used (pipes are
    used for the communication instead of stdin/stdout). Instead, there's a
    helper program `grale_socket_to_mpi.py`, for example:

    .. code-block:: bash

        grale_socket_to_mpi.py 9998 mpirun grale_massdens_mpi

    """

    def __init__(self, port, feedbackObject = None):
        """Connect to localhost (127.0.0.1) on the specified `port` number. The
        idea is to use SSH tunnels to connect to another machine."""
        super(NetcatLensPlaneRenderer, self).__init__(port, feedbackObject, "LENSPLANE")

_defaultMassRenderer = [ None ]
_defaultLensRenderer = [ None ]

def getDefaultLensPlaneRenderer():
    """Returns the default lens plane renderer that was set using :func:`setDefaultLensPlaneRenderer`."""
    return _defaultLensRenderer[0]

def getDefaultMassRenderer():
    """Returns the default mass renderer that was set using :func:`setDefaultMassRenderer`."""
    return _defaultMassRenderer[0]

def setDefaultLensPlaneRenderer(x):
    """Sets the default lens plane renderer, see the module documentation for allowed values."""
    _defaultLensRenderer[0] = x

def setDefaultMassRenderer(x):
    """Sets the default mass renderer, see the module documentation for allowed values."""
    _defaultMassRenderer[0] = x

