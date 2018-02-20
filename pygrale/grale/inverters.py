"""TODO:"""

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

real_print = print
def print(*args): # Do a flush afterwards
    real_print(*args)
    sys.stdout.flush()

if "INVERTER_USE_TIMEDIO" in os.environ:
    print("Using timed IO")
    from . import timedio as timed_or_untimed_io
else:
    #print("Using untimed IO")
    from . import untimedio as timed_or_untimed_io

class InverterException(Exception):
    """TODO:"""
    pass

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

    def invert(self, moduleName, populationSize, gaParams, gridLensInversionParameters, returnNds):

        errFile = tempfile.TemporaryFile("w+t")
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
            proc = subprocess.Popen(self.args, stdin = subprocess.PIPE, stdout = subprocess.PIPE, env = env)
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

            factoryParams = gridLensInversionParameters.toBytes()
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
                line = io.readLine(1000000)
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
            for i in range(numSols):
                line = io.readLine(30)
                if not line.startswith(resultStr):
                    raise InverterException("Unexpected identifier from inverter process: '{}'".format(line))

                numBytes = int(line[len(resultStr):])
                # Read the resulting lens
                lensData = io.readBytes(numBytes)

                line = io.readLine(30)
                if not line.startswith(fitId):
                    raise InverterException("Unexpected identifier from inverter process: '{}'".format(line))
                fitnessValues = [ float(x) for x in line[len(fitId):].strip().split(" ") ]

                sols.append((lenses.GravitationalLens.fromBytes(lensData), fitnessValues))

            io.writeLine("EXIT")
            #time.sleep(10)

            #print("Waiting...")
            proc.wait()
            #print("Done waiting")

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

def _usageAndDefaultParamsHelper(moduleName, exeName, key):
 
    try:
        proc = subprocess.Popen([exeName], stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    except Exception as e:
        raise InverterException("Unable to start process '{}': {}".format(exeName, e)) 

    inFd = proc.stdin.fileno()
    outFd = proc.stdout.fileno()
    io = timed_or_untimed_io.IO(outFd, inFd)

    line = io.readLine(30)
    invId = "GAINVERTER:"
    if not line.startswith(invId):
        raise InverterException("Unexpected idenficiation from process '{}': '{}'".format(exeName, line))

    version = line[len(invId):]

    io.writeLine("MODULE:" + moduleName)

    line = io.readLine(10)
    if not line.startswith(key + ":"):
        raise InverterException("Process '{}' didn't respond with expected '{}'".format(exeName, key))

    idx = line.find(":")
    numBytes = int(line[idx+1:])
    
    retVal = None
    if numBytes > 0:
        retVal = os.read(outFd, numBytes)

    io.writeLine("EXIT")

    return retVal

def getInversionModuleUsage(moduleName):
    """TODO:"""
    usageBytes = _usageAndDefaultParamsHelper(moduleName, "grale_invert_usage", "USAGE")
    return usageBytes.decode()

def getInversionModuleDefaultConfigurationParameters(moduleName):
    """TODO:"""
    confBytes = _usageAndDefaultParamsHelper(moduleName, "grale_invert_confparamdefaults", "CONFPARAMDEFAULTS")
    if confBytes:
        return inversionparams.ConfigurationParameters.fromBytes(confBytes).asDict()

    return None

def calculateFitness(moduleName, inputImages, zd, fitnessObjectParameters, lens):
    fitness, description = None, None

    try:
        proc = subprocess.Popen(["grale_invert_calcfitness"], stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    except Exception as e:
        raise InverterException("Unable to start process '{}': {}".format(exeName, e)) 

    inFd = proc.stdin.fileno()
    outFd = proc.stdout.fileno()
    io = timed_or_untimed_io.IO(outFd, inFd)

    line = io.readLine(30)
    invId = "GAINVERTER:"
    if not line.startswith(invId):
        raise InverterException("Unexpected idenficiation from process '{}': '{}'".format(exeName, line))

    version = line[len(invId):]

    io.writeLine("MODULE:" + moduleName)

    # Abusing the GridLensInversionParameters for this
    from . import grid
    g = grid.createUniformGrid(1, [0,0], 1) # Just a dummy grid
    g = grid._fractionalGridToRealGrid(g)
    Dd = 1.0 if not lens else lens.getLensDistance()
    params = inversionparams.GridLensInversionParameters(1, inputImages, g, Dd, zd, 1.0, baseLens = lens, 
                                                         fitnessObjectParameters=fitnessObjectParameters)

    factoryParams = params.toBytes()
    factoryParamsLen = len(factoryParams)
    io.writeLine("GAFACTORYPARAMS:{}".format(factoryParamsLen))
    io.writeBytes(factoryParams)

    while True:
        line = io.readLine(60*60*24*365*1000)
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

    io.writeLine("EXIT")

    return (fitness, description)

class SingleProcessInverter(Inverter):
    """TODO:"""
    def __init__(self, feedbackObject = None):
        """TODO:"""
        super(SingleProcessInverter, self).__init__([ "grale_invert_single" ], "Single process", feedbackObject=feedbackObject)

class MPIProcessInverter(Inverter):
    """TODO:"""
    def __init__(self, numProcesses = None, feedbackObject = None):
        """TODO:"""

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
    """TODO:"""
    def __init__(self, ipString, portNumber, feedbackObject = None):
        """TODO:"""
        super(ClientServerProcessInverter, self).__init__([ "grale_invert_clientserver", ipString, str(portNumber)],
                                                          "Client-server (IP {}, port {})".format(ipString, portNumber),
                                                          feedbackObject=feedbackObject)

class MPICSProcessInverter(Inverter):
    """TODO:"""
    def __init__(self, numProcesses = None, feedbackObject = None, serverDebugLevel = 0):
        """TODO:"""

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
        """TODO:"""
        if self.proc:
            try:
                privutil.terminateProcess(self.proc, feedbackObject = self.feedback)
            except Exception as e:
                print("Ignoring exception when terminating mpigaserver: " + str(e))

    def __del__(self):
        self.destroy()

class LocalCSProcessInverter(Inverter):
    """TODO:"""
    def __init__(self, numProcesses = None, feedbackObject = None, serverHelperDebugLevel = 0):
        """TODO:"""

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
        """TODO:"""
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
    return _defaultInverter[0]

def setDefaultInverter(x):
    _defaultInverter[0] = x

