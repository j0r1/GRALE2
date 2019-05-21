from __future__ import print_function
from . import renderers
import random
import subprocess
import time
import os
import tempfile
import errno

class RenderNameException(Exception):
    pass

def _matchRendererToFeedbackObject(renderer, feedbackObject):

    from . import feedback

    obj = None
    if feedbackObject is None:
        if renderer is not None:
            obj = renderer.getFeedbackObject()

    else:
        # feedbackObject is not None

        if type(feedbackObject) == str: # Description of which class to use
            feedbackClass = feedback.getFeedbackClass(feedbackObject)
            obj = feedbackClass()
        else:
            obj = feedbackObject

        if renderer is not None:
            # Make sure the renderer uses the specified feedback object
            renderer.setFeedbackObject(obj)

    if obj is None:
        obj = feedback.Feedback() # Just a dummy one, so that the functions can be called

    return obj

def _getMassRenderer(renderNameOrObject):
    if renderNameOrObject == "default":
        renderNameOrObject = renderers.getDefaultMassRenderer()

    if renderNameOrObject is None:
        return None

    if type(renderNameOrObject) == str:
        if renderNameOrObject.lower() == "threads":
            return renderers.ThreadsMassDensityRenderer()

        if renderNameOrObject.lower() == "mpi":
            return renderers.MPIMassDensityRenderer()

        raise RenderNameException("Supported renderer names are 'Threads' and 'MPI'")

    return renderNameOrObject

def _getLensPlaneRenderer(renderNameOrObject):

    if renderNameOrObject == "default":
        renderNameOrObject = renderers.getDefaultLensPlaneRenderer()

    if renderNameOrObject is None:
        return None

    if type(renderNameOrObject) == str:
        if renderNameOrObject.lower() == "threads":
            return renderers.ThreadsLensPlaneRenderer()

        if renderNameOrObject.lower() == "mpi":
            return renderers.MPILensPlaneRenderer()

        if renderNameOrObject.lower() == "opencl":
            return renderers.OpenCLLensPlaneRenderer()

        raise RenderNameException("Supported renderer names are 'Threads', 'OpenCL' and 'MPI'")

    return renderNameOrObject

def initRendererAndFeedback(renderer, feedbackObject, renderClass):
    if renderClass == "MASSDENS":
        renderer = _getMassRenderer(renderer)
    elif renderClass == "LENSPLANE":
        renderer = _getLensPlaneRenderer(renderer)
    else:
        raise RenderNameException("Invalid renderClass value '{}'".format(renderClass))

    feedbackObject = _matchRendererToFeedbackObject(renderer, feedbackObject)
    return (renderer, feedbackObject)

def generateRandomIdentifier():
    chars = "abcdefghijklmnopqrstuvwxyz"
    return ''.join([ chars[int(random.random()*len(chars))] for i in range(20) ])

def initInverterAndFeedback(inverter, feedbackObject):
    from . import inverters
    from . import feedback

    if inverter == "default":
        inverter = inverters.getDefaultInverter()
    
    if type(inverter) == str:
        if inverter.lower() == "singlecore":
            inverter = inverters.SingleProcessInverter()
        elif inverter.lower() == "gdb":
            inverter = inverters.SingleProcessGdbInverter()
        elif inverter.lower() == "mpi":
            inverter = inverters.MPIProcessInverter()
        elif inverter.lower() == "localcs":
            inverter = inverters.LocalCSProcessInverter()
        elif inverter.lower() == "mpics":
            inverter = inverters.MPICSProcessInverter()
        elif inverter.lower() == "clientserver":
            inverter = inverters.ClientServerProcessInverter()

    obj = None
    if feedbackObject is None:
        obj = inverter.getFeedbackObject()

    else:
        if type(feedbackObject) == str: # Description of which class to use
            feedbackClass = feedback.getFeedbackClass(feedbackObject)
            obj = feedbackClass()
        else:
            obj = feedbackObject
        
        inverter.setFeedbackObject(obj)

    if obj is None:
        obj = feedback.Feedback() # Just a dummy one, so that the functions can be called

    return (inverter, obj)

def initCosmology(cosm):
    from . import cosmology
    if cosm == "default":
        cosm = cosmology.getDefaultCosmology()
    
    return cosm

_terminateInfo = {
    "iswin": None,
    "havepkill": None,
    "havepgrep": None
}

if _terminateInfo["iswin"] is None:
    if hasattr(subprocess, 'STARTUPINFO'):
        _terminateInfo["iswin"] = True
        _terminateInfo["havepkill"] = False
        _terminateInfo["havepgrep"] = False
    else:
        _terminateInfo["iswin"] = False
        _terminateInfo["havepkill"] = False
        _terminateInfo["havepgrep"] = False
        if "PATH" in os.environ:
            for p in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(p, "pkill")):
                    _terminateInfo["havepkill"] = True
                if os.path.exists(os.path.join(p, "pgrep")):
                    _terminateInfo["havepgrep"] = True


def terminateProcess(proc, feedbackObject = None, maxHalfTime = 1.0):

    if feedbackObject is None:
        from . import feedback
        feedbackObject = feedback.Feedback()

    isWin = _terminateInfo["iswin"]
    havePKill = _terminateInfo["havepkill"]
    havePGrep = _terminateInfo["havepgrep"]

    subProcIds = [ ]

    if isWin:
        feedbackObject.onStatus("On windows, will use taskkill for process termination")
    if havePKill:
        feedbackObject.onStatus("Pkill is available for process termination")

    if havePGrep:
        subProcIds = list(map(int, subprocess.check_output([ "pgrep", "-P", str(proc.pid)]).decode().splitlines()))

    # First terminate
    if proc.poll() is None:
        feedbackObject.onStatus("Terminating process")
        t0 = time.time()
        while time.time() - t0 < maxHalfTime:
            try:
                if not isWin:
                    if havePKill:
                        subprocess.call(["pkill", "-P", str(proc.pid)])

                    proc.terminate()
                    time.sleep(0.05)
                else:
                    subprocess.call(["taskkill", "/F", "/T", "/PID", str(proc.pid)])
            except Exception as e:
                feedbackObject.onStatus("Error terminating background process: " + str(e))

            if proc.poll() is not None:
                break

            time.sleep(0.05)

    # Then kill
    if proc.poll() is None:
        feedbackObject.onStatus("Killing process")
        t0 = time.time()
        while time.time() - t0 < maxHalfTime:
            try:
                proc.kill()
                time.sleep(0.05)
            except Exception as e:
                feedbackObject.onStatus("Error killing background process: " + str(e))
                pass

            if proc.poll() is not None:
                break

            time.sleep(0.05)

    # Check child processes again

    if havePGrep:
        startTime = time.time()	
        while subProcIds and time.time()-startTime < 10.0:
            time.sleep(0.05)
            for opts in [ [], [ "-KILL" ] ]:
                time.sleep(0.1)
                newProcIds = [ ]
                for pid in subProcIds:
                    if _pid_exists(pid):
                        newProcIds.append(pid)

                        cmd = [ "kill" ] + opts + [ str(pid) ]
                        #print(cmd)
                        subprocess.call(cmd)

                subProcIds = newProcIds

# From https://stackoverflow.com/a/6940314/2828217
def _pid_exists(pid):
    """Check whether pid exists in the current process table.
    UNIX only.
    """
    if pid < 0:
        return False
    if pid == 0:
        # According to "man 2 kill" PID 0 refers to every process
        # in the process group of the calling process.
        # On certain systems 0 is a valid PID but we have no way
        # to know that in a portable fashion.
        raise ValueError('invalid PID 0')
    try:
        os.kill(pid, 0)
    except OSError as err:
        if err.errno == errno.ESRCH:
            # ESRCH == No such process
            return False
        elif err.errno == errno.EPERM:
            # EPERM clearly means there's a process to deny access to
            return True
        else:
            # According to "man 2 kill" possible error values are
            # (EINVAL, EPERM, ESRCH)
            raise
    else:
        return True

class PipePair(object):
    def __init__(self):

        self.rdFileName = tempfile.NamedTemporaryFile().name
        self.wrFileName = tempfile.NamedTemporaryFile().name
        os.mkfifo(self.rdFileName)
        os.mkfifo(self.wrFileName)
        self.rdPipeDesc = os.open(self.rdFileName, os.O_RDWR)
        self.wrPipeDesc = os.open(self.wrFileName, os.O_RDWR)

        #print("PIPE FILE NR IS:", self.rdPipeDesc, self.wrPipeDesc)

    def __del__(self):
        try:
            os.close(self.rdPipeDesc)
        except:
            pass
        try:
            os.unlink(self.rdFileName)
        except:
            pass
        try:
            os.close(self.wrPipeDesc)
        except:
            pass
        try:
            os.unlink(self.wrFileName)
        except:
            pass


def _wait(proc, maxtime):
    startTime = time.time()
    while time.time() - startTime < maxtime:
        if proc.poll() is not None:
            break
        time.sleep(0.1)

