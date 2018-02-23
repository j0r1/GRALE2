from __future__ import print_function
from . import renderers
import random
import subprocess
import time
import os
import tempfile

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
        if renderNameOrObject.lower() == "openmp":
            return renderers.OpenMPMassDensityRenderer()

        if renderNameOrObject.lower() == "mpi":
            return renderers.MPIMassDensityRenderer()

        raise RenderNameException("Supported renderer names are 'OpenMP' and 'MPI'")

    return renderNameOrObject

def _getLensPlaneRenderer(renderNameOrObject):

    if renderNameOrObject == "default":
        renderNameOrObject = renderers.getDefaultLensPlaneRenderer()

    if renderNameOrObject is None:
        return None

    if type(renderNameOrObject) == str:
        if renderNameOrObject.lower() == "openmp":
            return renderers.OpenMPLensPlaneRenderer()

        if renderNameOrObject.lower() == "mpi":
            return renderers.MPILensPlaneRenderer()

        if renderNameOrObject.lower() == "opencl":
            return renderers.OpenCLLensPlaneRenderer()

        raise RenderNameException("Supported renderer names are 'OpenCL', 'OpenMP' and 'MPI'")

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

_terminateInfo = {
    "iswin": None,
    "havepkill": None
}

if _terminateInfo["iswin"] is None:
    if hasattr(subprocess, 'STARTUPINFO'):
        _terminateInfo["iswin"] = True
        _terminateInfo["havepkill"] = False
    else:
        _terminateInfo["iswin"] = False
        _terminateInfo["havepkill"] = False
        if "PATH" in os.environ:
            for p in os.environ["PATH"].split(os.pathsep):
                if os.path.exists(os.path.join(p, "pkill")):
                    _terminateInfo["havepkill"] = True

def terminateProcess(proc, feedbackObject = None, maxHalfTime = 1.0):

    if feedbackObject is None:
        from . import feedback
        feedbackObject = feedback.Feedback()

    isWin = _terminateInfo["iswin"]
    havePKill = _terminateInfo["havepkill"]

    if isWin:
        feedbackObject.onStatus("On windows, will use taskkill for process termination")
    if havePKill:
        feedbackObject.onStatus("Pkill is available for process termination")

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
            except Exception as e:
                feedbackObject.onStatus("Error killing background process: " + str(e))
                pass

            if proc.poll() is not None:
                break

            time.sleep(0.05)


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

