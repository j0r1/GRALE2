"""This module provides different 'inverters' that can be used to execute
the genetic algorithm for lens inversion. Depending on the platform you're
working on, you may be able to use MPI to distribute the needed calculations
(typically speeding up the process) for example. 

In the :mod:`inversion` module, you can specify the inverter to use as
a (case insensitive) string, or as an instance of the available inverter
classes. Available are:

 - "Threads"/"Threads:numthreads" or an instance of :class:`ThreadsInverter`
 - "MPI"/"MPI:numprocesses"/"MPI:numprocesses:numthreadsperprocess:extraoptions" or an
   instance of :class:`MPIProcessInverter`

"""
from __future__ import print_function
from . import lenses
from . import inversionparams
from . import privutil
import subprocess
import os
import tempfile
import time
import struct
import socket
import sys

debugOutput = False
debugDirectStderr = False
debugCaptureProcessCommandsFile = None
communicationTimeout = 300

real_print = print
def print(*args): # Do a flush afterwards
    real_print(*args)
    sys.stdout.flush()

# For now, always use the timed version: better chance of detecting a
# crashed MPI process (since we're using a named pipe in that case, the
# crashed subprocess doesn't seem to register as a closed connection?
# TODO: alternatively the 'mpi' version could be disabled, and 'mpics'
#       used?
if not hasattr(subprocess, 'STARTUPINFO'): # Not Windows: can't use 'poll' there
    # print("Using timed IO")
    from . import timedio as timed_or_untimed_io
else:
    #print("Using untimed IO")
    from . import untimedio as timed_or_untimed_io

class InverterException(Exception):
    """An exception that's generated if something goes wrong in this module."""
    pass

class IOWrapper(timed_or_untimed_io.IO):
    def __init__(self, outFD, inFD, fname):
        self.f = open(fname, "wb")
        super(IOWrapper, self).__init__(outFD, inFD)

    def writeBytes(self, b):
        self.f.write(b)
        return super(IOWrapper, self).writeBytes(b)

def _writeParameters(io, paramName, paramBytes):
    io.writeLine("{}:{}".format(paramName,len(paramBytes)))
    io.writeBytes(paramBytes)

def _readString(io, prefix, operation = lambda x: x):
    prefix += ":"
    line = io.readLine(30)
    if not line.startswith(prefix):
        raise InverterException("Unexpected identifier from inverter process: '{}'".format(line))
    return operation(line[len(prefix):])

def _readInteger(io, prefix):
    return _readString(io, prefix, lambda x: int(x))

def _startProcess(args, extraEnv):
    try:
        errFile = tempfile.TemporaryFile("w+t") if not debugDirectStderr else None
        env = privutil._mergeExtraEnvironmentVariables(extraEnv)
        proc = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = subprocess.PIPE, env = env, stderr = errFile)
        return (proc, errFile)
    except Exception as e:
        raise InverterException("Unable to start inversion process: {}".format(e))

def _waitForInverterIdentification(io):
    line = io.readLine(30)
    invId = "GAINVERTER:"
    if not line.startswith(invId):
        raise InverterException("Unexpected idenficiation from inverter process: '{}'".format(line))

    version = line[len(invId):]
    return version

class Inverter(object):
    def __init__(self, args, inversionType, extraEnv = None, feedbackObject = None, readDescriptor = None, 
                 writeDescriptor = None):
        self.invType = inversionType
        self.args = args
        self.extraEnv = extraEnv
        self.feedback = feedbackObject
        self.readDescriptor = readDescriptor
        self.writeDescriptor = writeDescriptor

    def setDescriptors(self, readDescriptor, writeDescriptor):
        self.readDescriptor = readDescriptor
        self.writeDescriptor = writeDescriptor

    def setFeedbackObject(self, feedbackObject):
        self.feedback = feedbackObject

    def getFeedbackObject(self):
        return self.feedback

    def onStartedProcess(self, proc):
        pass

    def _getIO(self, proc):
        inFd = proc.stdin.fileno() if self.writeDescriptor is None else self.writeDescriptor
        outFd = proc.stdout.fileno() if self.readDescriptor is None else self.readDescriptor

        io = timed_or_untimed_io.IO(outFd, inFd) if not debugCaptureProcessCommandsFile else IOWrapper(outFd, inFd, debugCaptureProcessCommandsFile)
        return io

    def _writeInverterParameters(self, io, moduleName, calcType, populationSize,
                                 gaParams, lensInversionParameters, convergenceParams):
        io.writeLine("FITNESSOBJECT:{}".format(moduleName))
        io.writeLine("CALCULATOR:{}".format(calcType))
        io.writeLine("POPULATIONSIZE:{}".format(populationSize))

        if not gaParams: gaParams = { }
        gaParamsBytes = inversionparams.GAParameters(**gaParams).toBytes()

        _writeParameters(io, "GAPARAMS", gaParamsBytes)
        _writeParameters(io, "GAFACTORYPARAMS", lensInversionParameters.toBytes())
        _writeParameters(io, "CONVERGENCEPARAMS", convergenceParams.toBytes())

    def _run(self, io, returnNds):
        
        io.writeLine("RUN" if not returnNds else "RUN:NDS")
        
        # TODO: currently, status or progress are not used
        statusStr = "STATUS:"
        critId = "CRITERIA:"

        while True:
            line = io.readLine(communicationTimeout) # Wait at most five minutes
            #print(line)
            if line.startswith(statusStr):
                s = line[len(statusStr):]
                self.onStatus(s)
            elif line.startswith(critId):
                break
            else:
                print(line)
                #raise InverterException("Unexpected line '{}'".format(line))
            
        fitnessCriteria = line[len(critId):].strip().split(" ")
        return fitnessCriteria

    def invert(self, moduleName, calcType, populationSize, gaParams, lensInversionParameters, returnNds, convergenceParams):

        proc, errFile = _startProcess(self.args, self.extraEnv)

        try:
            self.onStartedProcess(proc)
            self.onStatus("Using inversion process for type: " + self.invType)
            
            io = self._getIO(proc)
            version = _waitForInverterIdentification(io)
            self.onStatus("Version info: " + version)

            self._writeInverterParameters(io, moduleName, calcType, populationSize,
                                          gaParams, lensInversionParameters, convergenceParams)

            # Actually run the inversion
            fitnessCriteria = self._run(io, returnNds)

            numSols = _readInteger(io, "NUMSOLS")
            self.onStatus("Expecting {} solutions".format(numSols))

            sols = []
            for _ in range(numSols):
                numBytes = _readInteger(io, "RESULT")
                lensData = io.readBytesUntimed(numBytes) # Read the resulting lens
                fitnessValues = _readString(io, "FITNESS", lambda l: [ float(x) for x in l.strip().split(" ") ])

                sols.append((lenses.GravitationalLens.fromBytes(lensData), fitnessValues))

            # print("Writing EXIT")
            io.writeLine("EXIT")

            self.onStatus("Waiting for process to finish...")
            privutil._wait(proc, 1)
        except InverterException:
            raise
        except Exception as e:
            raise InverterException(privutil._getErrorLineFromErrFile(errFile, debugOutput))
        finally:
            privutil._closeStdInOutAndterminateProcess(proc, self.feedback, debugOutput)
        
        if returnNds:
            return(sols, fitnessCriteria)

        if len(sols) != 1:
            raise InverterException("Internal error: expecting just one solution but got {}".format(len(sols)))

        return (sols[0][0], sols[0][1], fitnessCriteria)

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

def _commonModuleCommunication(moduleName, calcType, exeName, callback, **callbackArgs):

    proc, errFile = _startProcess([exeName], None)
    io = timed_or_untimed_io.IO(proc.stdout.fileno(), proc.stdin.fileno())

    try:
        _waitForInverterIdentification(io)

        io.writeLine("FITNESSOBJECT:" + moduleName)
        io.writeLine("CALCULATOR:" + calcType)

        retVal = callback(io, **callbackArgs)

        io.writeLine("EXIT")
    except Exception as e:
        raise InverterException(privutil._getErrorLineFromErrFile(errFile, debugOutput))

    return retVal

def _usageAndDefaultParamsHelper(moduleName, exeName, key):
 
    def f(io):
        line = io.readLine(10)
        if not line.startswith(key + ":"):
            raise InverterException("Process '{}' didn't respond with expected '{}'".format(exeName, key))

        numBytes = int(line[line.find(":")+1:])
        return None if numBytes <= 0 else io.readBytes(numBytes)

    retVal = _commonModuleCommunication(moduleName, "(notused)", exeName, f)

    return retVal

# TODO: should this be documented?
def getInversionModuleUsage(moduleName):
    usageBytes = _usageAndDefaultParamsHelper(moduleName, "grale_invert_usage_new", "USAGE")
    return usageBytes.decode()

# TODO: should this be documented?
def getInversionModuleDefaultConfigurationParameters(moduleName):
    confBytes = _usageAndDefaultParamsHelper(moduleName, "grale_invert_confparamdefaults_new", "CONFPARAMDEFAULTS")
    if confBytes:
        return inversionparams.ConfigurationParameters.fromBytes(confBytes).asDict()

    return None

# if lens is None and bpImages is None, only the fitness description is returned
# if lens is set, it is used to backproject the images
# if bpImages is set, they are assumed to be the backprojected images
def calculateFitness(moduleName, inputImages, zd, fitnessObjectParameters, lens = None, bpImages = None):

    if lens is None and bpImages is None:
        typeStr = "fitnessdescription"
    elif lens is None and bpImages is not None:
        typeStr = "precalculated"
        if len(inputImages) != len(bpImages):
            raise InverterException("Backprojected images list should have same length as images list")
    elif lens is not None and bpImages is None:
        typeStr = "lens"
    else:
        raise InverterException("Can't handle 'lens' and 'bpImages' both bein set, at least one should be 'None'")

    def f(io):

        fitness, description = None, None

        io.writeLine("TYPE:" + typeStr)

        # Abusing the LensInversionParametersSinglePlaneCPU for this
        from . import grid
        g = grid.createUniformGrid(1, [0,0], 1) # Just a dummy grid
        g = grid._fractionalGridToRealGrid(g)
        Dd = 1.0 if not lens else lens.getLensDistance()
        params = inversionparams.LensInversionParametersSinglePlaneCPU(inputImages, { "gridSquares": g },
                                                             Dd, zd, 1.0, baseLens = None, 
                                                             fitnessObjectParameters=fitnessObjectParameters)

        factoryParams = params.toBytes()
        factoryParamsLen = len(factoryParams)
        io.writeLine("GAFACTORYPARAMS:{}".format(factoryParamsLen))
        io.writeBytes(factoryParams) 

        if lens:
            lensBytes = lens.toBytes()
            io.writeLine("LENS:{}".format(len(lensBytes)))
            io.writeBytes(lensBytes)
        elif bpImages:
            for img in bpImages:
                imgBytes = img.toBytes()
                io.writeLine("IMGDATA:{}".format(len(imgBytes)))
                io.writeBytes(imgBytes)

        while True:
            # This can take a while
            # TODO send some keepalive message so that this timeout can be lower
            line = io.readLine(60*60*24) 
            if line == "DONE":
                break

            p = line.split(":")
            if len(p) != 2:
                raise InverterException("Expecting two parts for separator ':', but got line '{}'".format(line))
            if p[0] == "FITDESC":
                description = p[1]
            elif p[0] == "FITNESS":
                fitness = list(map(float, p[1].split(",")))
            else:
                raise InverterException("Expecting FITDESC or FITNESS, but got line '{}'".format(line))

        return fitness, description

    fitness, description = _commonModuleCommunication(moduleName, "(notused)", "grale_invert_calcfitness_new", f)

    return (fitness, description)

class ThreadsInverter(Inverter):

    def __init__(self, numThreads = 0, feedbackObject = None):

        numThreads = _getNumHelpers(numThreads)
        super(ThreadsInverter, self).__init__([ "grale_invert_newga" ],
                                                     "Thread based inverter",
                                                     extraEnv = { "GRALE_NUMTHREADS": str(numThreads) },
                                                     feedbackObject=feedbackObject)

def _createRandomTmpPath(length = 16):
    import string
    import numpy as np

    tmpDir = tempfile.gettempdir()
    chars = string.ascii_letters + string.digits
    return os.path.join(tmpDir, "".join(map(lambda x: chars[x], np.random.randint(0,len(chars),length))))

def _createBoundUnixSocket(socketPath):
    import socket

    s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) 
    s.bind(socketPath)
    s.listen(1)
    return s

class MPIProcessInverter(Inverter):

    def __init__(self, numProcesses = None, numThreadsPerProcess = 1, extraArgs = [], incomingConnectionTimeout = 10, feedbackObject = None):

        # Using stdin/stdout does not seem to work well for MPI (at least not openmpi),
        # so we'll use a UNIX socket        
        self.socketPath = _createRandomTmpPath()
        
        # We'll create the socket here, wait for an incoming connection in onStartedProcess
        self.servSock = _createBoundUnixSocket(self.socketPath)
        self.incomingConnectionTimeout = incomingConnectionTimeout

        npArgs = [ ] if numProcesses is None else [ "-np", str(numProcesses) ]

        if numThreadsPerProcess > 1:
            # TODO: use feedback object?
            print("""
WARNING: requesting more than one thread per MPI process. Depending on the
         MPI implementation it's possible that extra options need to be set 
         to be able to use the threads effectively. E.g. with openmpi the 
         option "--bind-to socket" or "--bind-to none" may be necessary.
         """)
            if len(extraArgs) == 0: # Try to detect if it's OpenMPI and add --bind-to socket automatically
                try:
                    out = subprocess.check_output(["mpirun", "--version"]).decode()
                    if "Open MPI" in out:
                        print("Detected OpenMPI, adding --bind-to socket")
                        extraArgs += [ "--bind-to", "socket" ]
                except Exception as e:
                    print("WARNING: unable to detect OpenMPI version/implementation:", e)
            time.sleep(2)

        npArgs += extraArgs

        super(MPIProcessInverter, self).__init__( [ "mpirun" ] + npArgs + [ "grale_invert_newgampi", self.socketPath, str(numThreadsPerProcess) ], 
                                                  "MPI inverter", feedbackObject=feedbackObject)

    def onStartedProcess(self, proc):

        import select

        if self.feedback:
            self.feedback.onStatus("Waiting for incoming connection on UNIX socket " + self.socketPath)

        # Wait for connection
        rlist, _, _ = select.select([self.servSock], [], [], self.incomingConnectionTimeout)
        if not rlist:
            raise InverterException("Timeout while waiting for connection from MPI helper process")

        self.sock, _ = self.servSock.accept()
        self.servSock.close()

        self.setDescriptors(self.sock.fileno(), self.sock.fileno())

    def __del__(self):
        try:
            self.sock.close()
        except Exception as e:
            pass

        try:
            os.unlink(self.socketPath)
        except Exception as e:
            pass

# TODO: allow multiple threads as well?
class SingleProcessGdbInverter(Inverter):
    """If this inverter is used, a single process, single core method is used
    which is started using the debugger `GDB <https://www.gnu.org/software/gdb/>`_.
    To be able to interact with GDB, it is started in an `xterm <https://en.wikipedia.org/wiki/Xterm>`_
    process. This is meant to make debugging this multi-process architecure 
    somewhat easier."""
    def __init__(self, feedbackObject = None):
        """Initializes an instance of this class; a specific :mod:`feedback <grale.feedback>`
        object can be specified for status updates."""
        self.pipePair = privutil.PipePair() # Keep it for the lifetime of this object
        pp = self.pipePair
        super(SingleProcessGdbInverter, self).__init__([ "xterm", "-e", 
                                                         "gdb grale_invert_newga -ex 'set args {} {}' ; echo sleeping 10 seconds; sleep 10".format(pp.wrFileName, pp.rdFileName),
                                                         ], 
                                                         "Single process GDB", feedbackObject=feedbackObject,
                                                         extraEnv = { "GRALE_NUMTHREADS": "1" },
                                                         readDescriptor=pp.rdPipeDesc, 
                                                         writeDescriptor=pp.wrPipeDesc)

def _getNumHelpers(n):
    if n is None or n < 1:
        import multiprocessing
        return multiprocessing.cpu_count()
    return n

_defaultInverter = [ "threads" ]

def getDefaultInverter():
    """Returns the default inverter to use when running the genetic algorithm
    for the lens inversion. Can be changed using :func:`setDefaultInverter`."""
    return _defaultInverter[0]

def setDefaultInverter(x):
    """Sets the default inverter to use when running the genetic algorithm."""
    _defaultInverter[0] = x

def createInverterFromString(inverter):
    """Creates on of the inverters based on the specified inverter name.

    Can be one of (case insensitive)
    - "Threads"/"Threads:numthreads"
    - "MPI"/"MPI:numprocesses"/"MPI:numprocesses:numthreadsperprocess"

    """

    threadsPrefix = "threads:"
    mpiPrefix = "mpi:"

    if inverter.lower() == "threads":
        return ThreadsInverter()

    if inverter.lower() == "mpi":
        return MPIProcessInverter()

    if inverter.lower() == "gdb":
        return SingleProcessGdbInverter()

    if inverter.lower().startswith(mpiPrefix):
        rest = inverter[len(mpiPrefix):]
        parts = rest.split(":")
        extraArgs = []
        if len(parts) == 1:
            numThreadsPerProcess = 1
            numNodes = int(parts[0])
        else:
            numNodes = int(parts[0]) if parts[0] else None
            numThreadsPerProcess = int(parts[1])
            if len(parts) == 2:
                pass
            elif len(parts) == 3:
                extraArgs = parts[2].split(" ")
            else:
                raise InverterException("Unknown syntax for mpi inverter")
        return MPIProcessInverter(numNodes, numThreadsPerProcess, extraArgs)

    if inverter.lower().startswith(threadsPrefix):
        numThreads = int(inverter[len(threadsPrefix):])
        return ThreadsInverter(numThreads)

    raise InverterException("The specified string '{}' cannot be interpreted as an inverter".format(inverter))
