"""This module is used to provide you with feedback during some calculations
that can take a while, like rendering the mass density of a gravitational lens
in a distributed way using MPI.

When a feedback mechanism is required by one of the functions in :mod:`grale.plotutil`
you'll normally just specify a name for the feedback mechanism, which can be one of

 - ``"none"``: in this case no feedback will be generated
 - ``"stdout"``: feedback will be written line per line to standard output
 - ``"notebook"``: this mechanism can be used if an IPython notebook is being used, in
   which case e.g. a progress bar will be shown.
 - ``"default"``: selects the feedback mechanism that has been set by the :func:`setDefaultFeedback`
   function.

For example, in an IPython notebook you could do something like::

    import grale.feedback as feedback
    
    feedback.setDefaultFeedback("notebook")

and various rendering calculations will use a nice progress bar and a status message box.

Instead of specifying a feedback mechanism as a string, you can also specify the feedback
object to be used. This could be helpful if you wanted to create your own feedback mechanism
for example. In that case, you'll need to define an object that contains ``onStatus`` and
``onProgress`` member functions, like in the :class:`Feedback` base class.
"""

from __future__ import print_function
try:
    from IPython.display import display
except:
    def display(*args, **kwargs):
        pass

class Feedback(object):
    """The base class for feedback objects, which doesn't do anything by itself. It
    does define the :func:`onStatus` and :func:`onProgress` member functions that
    are required for a feedback object."""

    def onStatus(self, s):
        """This function will be called if a status message ``s`` needs to be displayed."""
        pass

    def onProgress(self, x):
        """This function will be called if a progress percentage ``x`` is reported."""
        pass

class StdoutFeedback(Feedback):
    """Derives from :class:`Feedback` and writes the status messages and percentages to
    the standard output."""

    def __init__(self):
        self.lastProgress = -1

    def onStatus(self, s):
        print("STATUS:", s)

    def onProgress(self, x):
        
        ix = int(round(x))
        if ix == self.lastProgress:
            return

        self.lastProgress = ix

        print("PROGRESS:", x)

class NotebookFeedback(Feedback):
    """Derives from :class:`Feedback` and writes the status messages and percentages to
    widgets that can be displayed in an IPython notebook."""

    def __init__(self, minValue = 0, maxValue = 100):

        import ipywidgets

        self.output = ipywidgets.Output()
        with self.output:
            self.progress = ipywidgets.widget_float.FloatProgress(min=minValue, max=maxValue)
            self.status = ipywidgets.widget_string.Text()

        self.displayedStatus = False
        self.displayedProgress = False

    def __del__(self):
        self.progress.close()
        self.status.close()
        self.output.clear_output()

    def onStatus(self, s):
        if not self.displayedStatus:
            self.displayedStatus = True
            display(self.status)

        self.status.value = "Status: {}".format(s)

    def onProgress(self, x):
        if not self.displayedProgress:
            self.displayedProgress = True
            display(self.progress)

        self.progress.value = x


defaultFeedback = { }
defaultFeedback["class"] = StdoutFeedback

class FeedbackException(Exception):
    """This exception is thrown in case of an error in this module."""
    pass

def getFeedbackClass(feedbackName):
    """
    Returns the class that corresponds to ``feedbackName``, which can be one of
     - ``"none"``: in this case :class:`Feedback` will be returned, so no output will be generated
     - ``"stdout"``: the :class:`StdoutFeedback` class will be returned
     - ``"notebook"``: the :class:`NotebookFeedback` class will be returned
     - ``"default"``: returns the class that was previously set by calling :func:`setDefaultFeedback`
    """

    if feedbackName.lower() == "default":
        return defaultFeedback["class"]

    if feedbackName.lower() == "stdout":
        return StdoutFeedback

    if feedbackName.lower() == "notebook":
        return NotebookFeedback

    if feedbackName.lower() == "none":
        return Feedback

    raise FeedbackException("Invalid feedback mechanism name, should be 'default', 'stdio' or 'notebook', but is '{}'".format(feedbackName))

def setDefaultFeedback(feedbackNameOrClass):
    """
    Sets the feedback mechamism that's used when specifying ``"default"``. The
    ``feedbackNameOrClass`` parameter can either be one of the strings that's recognized
    by the :func:`getFeedbackClass` function, or it can be the name of a class derived
    from :class:`Feedback`.
    """
    if type(feedbackNameOrClass) == str:
        defaultFeedback["class"] = getFeedbackClass(feedbackNameOrClass)
    else:
        defaultFeedback["class"] = feedbackNameOrClass


