"""This module provides different 'inverters' that can be used to execute
the genetic algorithm for lens inversion. Depending on the platform you're
working on, you may be able to use MPI to distribute the needed calculations
(typically speeding up the process) for example. 

In the :mod:`inversion` module, you can specify the inverter to use as
a (case insensitive) string, or as an instance of the available inverter
classes. Available are:

 - "SingleCore" or an instance of :class:`SingleProcessInverter`
 - "GDB" or and instance of :class:`SingleProcessGdbInverter` (for debugging code)
 - "MPI"/"MPI:numprocesses" or an instance of :class:`MPIProcessInverter`
 - "ClientServer" or an instance of :class:`ClientServerProcessInverter`
 - "LocalCS"/"LocalCS:numprocesses" or an instance of :class:`LocalCSProcessInverter`
 - "MPICS"/"MPICS:numprocesses" or an instance of :class:`MPICSProcessInverter`

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

def _getErrorLineFromErrFile(errFile):
    if not errFile:
        return "No error file was opened, can't get more specific error line"

    errFile.flush()
    errFile.seek(0)
    errData = errFile.read()
    if debugOutput:
        print("DEBUG: full contents of error log:")
        print(errData)

    errLines = [ l.strip() for l in errData.splitlines() if len(l.strip()) > 0 ]
    if len(errLines) > 0:
        errLine = "Last line from error log: " + errLines[-1]
    else:
        errLine = "Something went wrong, but don't know what"
    return errLine

class Inverter(object):
    def __init__(self, args, inversionType, extraEnv = None, feedbackObject = None, readDescriptor = None, 
                 writeDescriptor = None):
        self.invType = inversionType
        self.args = args
        self.extraEnv = extraEnv
        self.version = ""
        self.feedback = feedbackObject
        self.readDescriptor = readDescriptor
        self.writeDescriptor = writeDescriptor

    def setFeedbackObject(self, feedbackObject):
        self.feedback = feedbackObject

    def getFeedbackObject(self):
        return self.feedback

    def invert(self, moduleName, populationSize, gaParams, lensInversionParameters, returnNds):

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
            raise InverterException("Unable to start inversion process: {}".format(e)) 

        try:
            self.onStatus("Using inversion process for type: " + self.invType)

            inFd = proc.stdin.fileno() if self.writeDescriptor is None else self.writeDescriptor
            outFd = proc.stdout.fileno() if self.readDescriptor is None else self.readDescriptor
            io = timed_or_untimed_io.IO(outFd, inFd)

            line = io.readLine(30)
            invId = "GAINVERTER:"
            if not line.startswith(invId):
                raise InverterException("Unexpected idenficiation from inverter process: '{}'".format(line))

            self.version = line[len(invId):]
            self.onStatus("Version info: " + self.version)

            io.writeLine("MODULE:{}".format(moduleName))
            io.writeLine("POPULATIONSIZE:{}".format(populationSize))

            if not gaParams:
                gaParams = { }
            gaParamsBytes = inversionparams.GAParameters(**gaParams).toBytes()

            gaParamsSize = len(gaParamsBytes)
            io.writeLine("GAPARAMS:{}".format(gaParamsSize))
            io.writeBytes(gaParamsBytes)

            factoryParams = lensInversionParameters.toBytes()
            factoryParamsLen = len(factoryParams)
            io.writeLine("GAFACTORYPARAMS:{}".format(factoryParamsLen))
            io.writeBytes(factoryParams)

            if not returnNds:
                io.writeLine("RUN")
            else:
                io.writeLine("RUN:NDS")
            
            # TODO: currently, status or progress are not used
            statusStr = "STATUS:"
            resultStr = "RESULT:"
            critId = "CRITERIA:"
            numSolsStr = "NUMSOLS:"
            fitId = "FITNESS:"

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

            line = io.readLine(30)
            if not line.startswith(numSolsStr):
                raise InverterException("Unexpected identifier from inverter process: '{}'".format(line))
            numSols = int(line[len(numSolsStr):])

            sols = []
            # TODO status message?
            print("Expecting {} solutions from line '{}'".format(numSols, line.strip()))
            for i in range(numSols):
                line = io.readLine(30)
                print("Read line '{}'".format(line))
                if not line.startswith(resultStr):
                    raise InverterException("Unexpected identifier from inverter process: '{}'".format(line))

                numBytes = int(line[len(resultStr):])
                # Read the resulting lens
                lensData = io.readBytes(numBytes)

                line = io.readLine(30)
                print("Read line '{}'".format(line))
                if not line.startswith(fitId):
                    raise InverterException("Unexpected identifier from inverter process: '{}'".format(line))
                fitnessValues = [ float(x) for x in line[len(fitId):].strip().split(" ") ]

                sols.append((lenses.GravitationalLens.fromBytes(lensData), fitnessValues))

            print("Writing EXIT")
            io.writeLine("EXIT")

            print("Waiting for process to finish...")
            #proc.wait()
            privutil._wait(proc, 1)
            print("Done waiting")
        except InverterException:
            raise
        except Exception as e:
            errLine = _getErrorLineFromErrFile(errFile)
            raise InverterException(errLine)
        finally:
            try:
                proc.stdin.close()
                proc.stdout.close()
            except:
                pass
            time.sleep(0.1)
            try:
                privutil.terminateProcess(proc, feedbackObject = self.feedback)
            except Exception as e:
                if debugOutput:
                    print("Ignoring exception when terminating program: " + str(e))
        
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

def _commonModuleCommunication(moduleName, exeName, callback, **callbackArgs):

    errFile = tempfile.TemporaryFile("w+t") if not debugDirectStderr else None
    try:
        proc = subprocess.Popen([exeName], stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr=errFile)
    except Exception as e:
        raise InverterException("Unable to start process '{}': {}".format(exeName, e)) 

    inFd = proc.stdin.fileno()
    outFd = proc.stdout.fileno()
    io = timed_or_untimed_io.IO(outFd, inFd)

    try:
        line = io.readLine(30)
        invId = "GAINVERTER:"
        if not line.startswith(invId):
            raise InverterException("Unexpected idenficiation from process '{}': '{}'".format(exeName, line))

        version = line[len(invId):]

        io.writeLine("MODULE:" + moduleName)

        retVal = callback(io, **callbackArgs)

        io.writeLine("EXIT")
    except Exception as e:
        errLine = _getErrorLineFromErrFile(errFile)
        raise InverterException(errLine)

    return retVal

def _usageAndDefaultParamsHelper(moduleName, exeName, key):
 
    def f(io):
        line = io.readLine(10)
        if not line.startswith(key + ":"):
            raise InverterException("Process '{}' didn't respond with expected '{}'".format(exeName, key))

        idx = line.find(":")
        numBytes = int(line[idx+1:])
        
        retVal = None
        if numBytes > 0:
            retVal = io.readBytes(numBytes)

        return retVal

    retVal = _commonModuleCommunication(moduleName, exeName, f)

    return retVal

# TODO: this uses a more low level moduleName, should this be documented?
def getInversionModuleUsage(moduleName):
    usageBytes = _usageAndDefaultParamsHelper(moduleName, "grale_invert_usage", "USAGE")
    return usageBytes.decode()

# TODO: this uses a more low level moduleName, should this be documented?
def getInversionModuleDefaultConfigurationParameters(moduleName):
    confBytes = _usageAndDefaultParamsHelper(moduleName, "grale_invert_confparamdefaults", "CONFPARAMDEFAULTS")
    if confBytes:
        return inversionparams.ConfigurationParameters.fromBytes(confBytes).asDict()

    return None

# TODO: this uses a more low level moduleName, should this be documented?
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
        params = inversionparams.LensInversionParametersSinglePlaneCPU(1, inputImages, { "gridSquares": g },
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

    fitness, description = _commonModuleCommunication(moduleName, "grale_invert_calcfitness", f)

    return (fitness, description)

class SingleProcessInverter(Inverter):
    """If this inverter is used, a single process, single core method is used. For
    a very simple inversion problem this may still be the most performant though.
    It is always available, irrespective of the platform you're working on."""

    def __init__(self, feedbackObject = None):
        """Initializes an instance of this class; a specific :mod:`feedback <grale.feedback>`
        object can be specified for status updates."""

        super(SingleProcessInverter, self).__init__([ "grale_invert_single" ], "Single process", feedbackObject=feedbackObject)

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
                                                         "gdb grale_invert_single -ex 'set args {} {}' ; echo sleeping 10 seconds; sleep 10".format(pp.wrFileName, pp.rdFileName),
                                                         ], 
                                                         "Single process GDB", feedbackObject=feedbackObject, 
                                                         readDescriptor=pp.rdPipeDesc, 
                                                         writeDescriptor=pp.wrPipeDesc)

class MPIProcessInverter(Inverter):
    """If MPI is available for your platform, this inverter will use the MPI system
    to distribute the calculations over the available processes."""

    def __init__(self, numProcesses = None, feedbackObject = None):
        """Initializes an instance of this class, optionally setting the number of
        MPI processes to use explicitly to `numProcesses`. A specific :mod:`feedback <grale.feedback>`
        object can be specified for status updates."""

        # Using stdin/stdout does not seem to work well for MPI (at least not openmpi),
        # so we'll use a different set of pipes
        self.pipePair = privutil.PipePair() # Keep it for the lifetime of this object
        pp = self.pipePair

        npArgs = [ ] if numProcesses is None else [ "-np", str(numProcesses) ]
        super(MPIProcessInverter, self).__init__( [ "mpirun" ] + npArgs + [ "grale_invert_mpi", pp.wrFileName, pp.rdFileName ], 
                                                  "MPI inverter", feedbackObject=feedbackObject, 
                                                  readDescriptor = pp.rdPipeDesc, writeDescriptor = pp.wrPipeDesc)

def _getNumHelpers(n):
    if n is None or n < 1:
        import multiprocessing
        return multiprocessing.cpu_count()
    return n

class ClientServerProcessInverter(Inverter):
    """This inverter will connect to a ``gaserver`` instance from the
    `MOGAL <http://research.edm.uhasselt.be/jori/mogal/>`_ library, which should already
    have several ``gahelper`` instances connected to it. The calculations
    will be divided over the available helpers."""

    def __init__(self, ipString = "127.0.0.1", portNumber = 9999, feedbackObject = None):
        """Initializes an instance of this class, where you must specify the
        IP address and port number at which the ``gaserver`` process can be
        reached.  A specific :mod:`feedback <grale.feedback>`
        object can be specified for status updates."""

        super(ClientServerProcessInverter, self).__init__([ "grale_invert_clientserver", ipString, str(portNumber)],
                                                          "Client-server (IP {}, port {})".format(ipString, portNumber),
                                                          feedbackObject=feedbackObject)

class MPICSProcessInverter(Inverter):
    """This inverter will start an ``mpigaserver`` instance from
    the `MOGAL <http://research.edm.uhasselt.be/jori/mogal/>`_ library, and
    will connect to it. In turn the ``mpigaserver`` program will distribute 
    the calculations using MPI."""

    def __init__(self, numProcesses = None, feedbackObject = None, serverDebugLevel = 0):
        """Initializes an instance of this class, optionally setting the number of
        MPI processes to use explicitly to `numProcesses`. A specific :mod:`feedback <grale.feedback>`
        object can be specified for status updates. The verbosity of the output of the
        ``mpigaserver`` process can be controller using `serverDebugLevel`"""

        self.proc = None

        # Obtain a port number to use, let's hope it will still be valid in a few seconds
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(('0.0.0.0', 0))
        bindIp, serverPort = s.getsockname()
        del s
    
        super(MPICSProcessInverter, self).__init__([ "grale_invert_clientserver", "127.0.0.1", str(serverPort) ],
                                                   "MPI client-server", feedbackObject=feedbackObject)

        # Get the module directory
        from . import inversion
        n = inversion._getModuleName("general")
        moduleDir = inversion._getModuleDirectory(n)

        npArgs = [ ] if numProcesses is None else [ "-np", str(numProcesses) ]

        self.proc = subprocess.Popen(["mpirun" ] + npArgs + [ "mpigaserver", str(serverDebugLevel), str(serverPort), moduleDir ])
        time.sleep(2) # wait a short while before allowing client to connect

    def destroy(self):
        """Stops the ``mpigaserver`` program."""
        if self.proc:
            try:
                privutil.terminateProcess(self.proc, feedbackObject = self.feedback)
            except Exception as e:
                print("Ignoring exception when terminating mpigaserver: " + str(e))

    def __del__(self):
        self.destroy()

class LocalCSProcessInverter(Inverter):
    """This inverter uses a similar approach as :class:`ClientServerProcessInverter`,
    but also starts the ``gaserver`` process and a number of ``gahelper`` processes
    itself, on the same host. This can also help to speed up calculations, by dividing
    them over a number of processes on the same computer. It will not have the same
    performance as MPI, but may still be better than using the 
    :class:`single process <SingleProcessInverter>` approach, especially when the
    calculations become more involved."""

    def __init__(self, numProcesses = None, feedbackObject = None, serverHelperDebugLevel = 0):
        """Initializes an instance of this class, optionally setting the number of
        helper processes to use explicitly to `numProcesses`. A specific :mod:`feedback <grale.feedback>`
        object can be specified for status updates. The verbosity of the output of the
        ``gaserver`` process can be controller using `serverDebugLevel`"""

        self.procs = [ ] # Do this eary, in case __del__ is called sooner than expected

        numHelpers = _getNumHelpers(numProcesses)
        if numHelpers < 1 or numHelpers > 256:
            raise InverterException("The number of helper processes should be at least one, and at most 256")

        # Obtain a port number to use, let's hope it will still be valid in a few seconds
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(('0.0.0.0', 0))
        bindIp, serverPort = s.getsockname()
        del s

        super(LocalCSProcessInverter, self).__init__([ "grale_invert_clientserver", "127.0.0.1", str(serverPort) ], 
                                                      "Local client-server ({} helpers)".format(numHelpers),
                                                      feedbackObject=feedbackObject)


        from . import inversion
        n = inversion._getModuleName("general")
        moduleDir = inversion._getModuleDirectory(n)
        
        p = subprocess.Popen(["gaserver", str(serverHelperDebugLevel), str(serverPort), moduleDir ])
        self.procs.append(p)
        time.sleep(2) # wait a short while before starting to connect the helpers

        for i in range(numHelpers):
            p = subprocess.Popen(["gahelper", str(serverHelperDebugLevel), "127.0.0.1", str(serverPort), moduleDir ])
            self.procs.append(p)

        time.sleep(0.5) # Wait a short while so that all helpers are completely detected
    
    def destroy(self):
        """Stops the ``gaserver`` and ``gahelper`` programs."""
        for proc in self.procs:
            try:
                privutil.terminateProcess(proc, feedbackObject = self.feedback)
            except Exception as e:
                print("Ignoring exception when terminating gaserver or gahelper: " + str(e))

        self.procs = []

    def __del__(self):
        self.destroy()

_defaultInverter = [ "singlecore" ]

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
    - "SingleCore"
    - "GDB"
    - "MPI"/"MPI:numprocesses"
    - "ClientServer"
    - "LocalCS"/"LocalCS:numprocesses"
    - "MPICS"/"MPICS:numprocesses"

    """

    localCsNodesPrefix = "localcs:"
    mpiNodesPrefix = "mpi:"
    csMpiNodesPrefix = "mpics:"

    if inverter.lower() == "singlecore":
        return SingleProcessInverter()

    if inverter.lower() == "gdb":
        return SingleProcessGdbInverter()

    if inverter.lower() == "mpi":
        return MPIProcessInverter()

    if inverter.lower().startswith(mpiNodesPrefix):
        numNodes = int(inverter[len(mpiNodesPrefix):])
        return MPIProcessInverter(numNodes)

    if inverter.lower() == "localcs":
        return LocalCSProcessInverter()

    if inverter.lower().startswith(localCsNodesPrefix):
        numNodes = int(inverter[len(localCsNodesPrefix):])
        return LocalCSProcessInverter(numNodes)

    if inverter.lower() == "mpics":
        return MPICSProcessInverter()

    if inverter.lower().startswith(csMpiNodesPrefix):
        numNodes = int(inverter[len(csMpiNodesPrefix):])
        return MPICSProcessInverter(numNodes)

    if inverter.lower() == "clientserver":
        return ClientServerProcessInverter()

    raise InverterException("The specified string '{}' cannot be interpreted as an inverter".format(inverter))
