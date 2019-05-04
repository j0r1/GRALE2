"""This module defines various plotting utilities, as well as 
classes to help you make animations.
"""

from __future__ import print_function
from . import privutil
from . import privutilcython
from . import images
from . import feedback
from . import gridfunction
from . import multiplane
from .bytestring import B, S
from .constants import *
import numpy as np
import subprocess
import os
import tempfile
import time
import json
import copy

class LensInfoException(Exception):
    """An exception that's raised in case of an error with a LensInfo object"""
    pass

class DensInfo(object):
    """In case you'd like to use the ``plotDensity`` functions to create
    a plot based on a 2D NumPy array, you can pass an instance of this class
    to the plot function. This can come in handy when creating the difference,
    average or standard deviation of several mass maps from :class:`LensInfo`
    objects."""
    def __init__(self, densitypointsOrPixels, bottomleft, topright, isPixels = False):
        """Initialize the object so that the 2D NumPy array `densitypointsOrPixels`
        describes an area from `bottomleft` to `topright`. The `isPixels` flag
        informs the constructor if the specified data contains values at regularly
        spaced points or values inside pixels. As the pixels can be calculated
        from the points, but not the other way around, use points if you have
        the option."""
        self.d = { }
        if not isPixels:
            self.d["densitypoints"] = densitypointsOrPixels
            self.d["numy"] = densitypointsOrPixels.shape[0]-1
            self.d["numx"] = densitypointsOrPixels.shape[1]-1
        else:
            self.d["densitypixels"] = densitypointsOrPixels
            self.d["numy"] = densitypointsOrPixels.shape[0]
            self.d["numx"] = densitypointsOrPixels.shape[1]

        self.d["bottomleft"] = copy.copy(bottomleft)
        self.d["topright"] = copy.copy(topright)

    def getBottomLeft(self):
        """Returns the `bottomleft` argument that was passed in the constructor,
        as a NumPy array."""
        return np.array(self.d["bottomleft"], dtype=np.double)

    def getTopRight(self):
        """Returns the `topright` argument that was passed in the constructor,
        as a NumPy array."""
        return np.array(self.d["topright"], dtype=np.double)

    def getArea(self):
        """Returns a dictionary with ``topright`` and ``bottomleft`` entries
        as specified in the constructor. Can come in handy when passing these
        arguments to another constructor to describe the same area."""
        return { "topright" : self.getTopRight(), "bottomleft": self.getBottomLeft() }

    def getDensityPoints(self, renderer = "default", feedbackObject = "default"):
        """Returns the `densitypoints` NumPy array that was set in the constructor."""
        return self.d["densitypoints"]

    def getDensityPixels(self, renderer = "default", feedbackObject = "default"):
        """Returns the density pixels, based on the `densitypoints` from the
        constructor. The `densitypoints` describe the values *at* specific points
        laid out on a grid; the pixels convert these values so that each entry
        becomes the integrated value over a pixel.
        
        In an instance of the :class:`DensInfo` type, the `renderer` and `feedbackObject`
        arguments are not used; in an instance of the derived class :class:`LensInfo`,
        they can be used to specify the :ref:`renderer <renderers>` to use
        to speed up the calculation, and with `feedbackObject` you can specify a
        particular :ref:`feedback mechanism <feedback>`."""

        if "densitypixels" in self.d:
            return self.d["densitypixels"]

        massMap = self.getDensityPoints(renderer, feedbackObject)
        numX, numY = self.d["numx"]+1, self.d["numy"]+1

        pixels = 0.25*massMap[0:numY-1,0:numX-1] + 0.25*massMap[0:numY-1,1:numX] + 0.25*massMap[1:numY,0:numX-1] + 0.25*massMap[1:numY,1:numX]
        self.d["densitypixels"] = pixels
        return pixels

    def getXPointCoordinates(self, repeat=True):
        """Returns the X coordinates for the `densitypoints` grid, based on the
        `bottomleft` and `topright` values from the constructor."""
        tmp = np.linspace(self.d["bottomleft"][0], self.d["topright"][0], self.d["numx"]+1)
        if not repeat:
            return tmp
        tmp = tmp.reshape((1,-1))
        return np.repeat(tmp, self.d["numy"]+1, 0)

    def getYPointCoordinates(self, repeat=True):
        """Returns the Y coordinates for the `densitypoints` grid, based on the
        `bottomleft` and `topright` values from the constructor."""
        tmp = np.linspace(self.d["bottomleft"][1], self.d["topright"][1], self.d["numy"]+1)
        if not repeat:
            return tmp
        tmp = tmp.reshape((-1, 1))
        return np.repeat(tmp, self.d["numx"]+1, 1)

    def getNumXPoints(self):
        """Returns the number of points in the X direction"""
        return self.d["numx"]+1

    def getNumYPoints(self):
        """Returns the number of points in the Y direction"""
        return self.d["numy"]+1

    def getNumXPixels(self):
        """Returns the number of pixels in the X direction"""
        return self.d["numx"]

    def getNumYPixels(self):
        """Returns the number of pixels in the Y direction"""
        return self.d["numy"]

class LensInfo(DensInfo):
    """This type of object can be used in the ``plotImagePlane`` and ``plotDensity``
    functions, and is based on a single or multi lensplane scenario.
    """

    def __init__(self, lens = None, bottomleft = None, topright = None, center = [0, 0], size = None, numx = None, numy = None, numxy = 511,
                 Ds = None, Dds = None, zd = None, zs = None, cosmology = "default"):
        """Constructs an object of this type.

        Arguments:
         * `lens`: for a single lens plane setting, this should be an instance of
           a :class:`GravitationalLens <grale.lenses.GravitationalLens>` derived
           class. In case a multi-lensplane scenario is to be used, this should be
           an array of ``(lens, redshift)`` tuples.
         * `bottomleft` and `topright`, or `size` and `center`: using these arguments
           you can describe the area for which plot information should be calculated.
         * `numx` and `numy`, or `numxy`: the number of columns and rows (or both) of 
           the grid for which properties (deflection angles, densities, ...) should be 
           calculated.
         * `Ds` and `Dds`, or `zs` (and possibly `zd`): these are relevant for
           image plane calculations, for which the position of the source needs to
           be known. For a single lens plane setting, you can either set both
           `Ds` and `Dds`, or both `zs` and `zd`. In case a multi lens plane setting
           is used, you only need to specify `zs`. Note that when using the redshift
           values, a cosmological model must be known as well.
         * `cosmology`: the cosmological model to be used. This is necessary when
           redshifts are to be used.
        """
        self.d = { }
        if lens: 
            self.d["lens"] = copy.copy(lens) # Note: lens can be both an array of (lens,z) tuples, or just a lens

        if bottomleft is not None and topright is not None:
            if size is not None and center is not None:
                raise LensInfoException("Both bottomleft/topright and center/size are set, can't determine what area to use")

            self.d["bottomleft"] = copy.copy(bottomleft)
            self.d["topright"] = copy.copy(topright)
        elif size is not None and center is not None:
            self.d["bottomleft"] = [ center[0] - size/2.0, center[1] - size/2.0 ]
            self.d["topright"] = [ center[0] + size/2.0, center[1] + size/2.0 ]
        else:
            raise LensInfoException("Either bottomleft/topright or center/size must be set")

        if numxy:
            self.d["numx"] = numxy
            self.d["numy"] = numxy

        if numx is not None: self.d["numx"] = numx
        if numy is not None: self.d["numy"] = numy

        if Ds is not None: self.d["Ds"] = Ds
        if Dds is not None: self.d["Dds"] = Dds
        if zd is not None: self.d["zd"] = zd
        if zs is not None: self.d["zs"] = zs

        cosmology = privutil.initCosmology(cosmology)
        self.cosm = copy.copy(cosmology)

    @staticmethod
    def fromLensPlane(lensPlane):
        """Creates a LensInfo instance based on a previously created :class:`LensInfo <grale.images.LensInfo>]`
        instance."""
        lens = lensPlane.getLens()
        ri = lensPlane.getRenderInfo()
        bl, tr = ri["bottomleft"], ri["topright"]
        xpoints, ypoints = ri["xpoints"], ri["ypoints"]
        li = LensInfo(lens, bottomleft = bl, topright = tr, numx = xpoints-1, numy = ypoints-1)
        li.d["lensplane"] = lensPlane
        return li

    def setSourceRedshift(self, zs):
        """Sets source redshift to `zs`."""
        self.d["zs"] = zs
        if "imageplane" in self.d:
            del self.d["imageplane"]

    def setSourceDistances(self, Ds, Dds):
        """Sets the source distances to `Ds` and `Dds`."""
        self.d["Ds"] = Ds
        self.d["Dds"] = Dds
        if "imageplane" in self.d:
            del self.d["imageplane"]

    def getLens(self):
        """Returns the lens (or multiple lenses) that was set at construction
        time."""
        if not "lens" in self.d:
            raise LensInfoException("No lens has been set")
        return self.d["lens"]

    def getLensPlane(self, renderer = "default", feedbackObject = "default"):
        """Returns the :class:`LensPlane` or :class:`MultiLensPlane` instance that
        corresponds to the settings from the constructor, performing the necessary
        calculations if needed. With the `renderer` argument, you can specify a specific :ref:`renderer <renderers>`
        to speed up the calculation, and with `feedbackObject` you can specify a
        particular :ref:`feedback mechanism <feedback>`."""

        renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "LENSPLANE")

        if "lensplane" in self.d:
            feedbackObject.onStatus("Reusing lens plane")
            return self.d["lensplane"]

        if not "lens" in self.d:
            raise LensInfoException("No lens available to calculate lens plane mapping from")

        if type(self.d["lens"]) != list:
            lensPlane = images.LensPlane(self.d["lens"], self.d["bottomleft"], self.d["topright"], self.d["numx"]+1, self.d["numy"]+1, renderer = renderer, feedbackObject = feedbackObject)
        else: # multiple lens
            if self.cosm is None:
                raise LensInfoException("For a multiple lensplane scenario, the cosmology must be set")

            lensPlane = multiplane.MultiLensPlane(self.d["lens"], self.d["bottomleft"], self.d["topright"], self.d["numx"]+1, self.d["numy"]+1, renderer = renderer, feedbackObject = feedbackObject, cosmology = self.cosm)

        self.d["lensplane"] = lensPlane
        return lensPlane

    def getImagePlane(self, renderer = "default", feedbackObject = "default"):
        """Returns the :class:`ImagePlane` or :class:`MultiImagePlane` instance that
        corresponds to the settings from the constructor, performing the necessary
        calculations if needed. With the `renderer` argument, you can specify a specific :ref:`renderer <renderers>`
        to speed up the calculation, and with `feedbackObject` you can specify a
        particular :ref:`feedback mechanism <feedback>`."""

        renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "LENSPLANE")

        if "imageplane" in self.d:
            feedbackObject.onStatus("Reusing image plane")
            return self.d["imageplane"]

        lensPlane = self.getLensPlane(renderer, feedbackObject)

        if type(self.d["lens"]) != list:
            if "Ds" in self.d and "Dds" in self.d:
                Ds = self.d["Ds"]
                Dds = self.d["Dds"]
            elif "zd" in self.d and "zs" in self.d and self.cosm:
                Dd = self.cosm.getAngularDiameterDistance(self.d["zd"])
                Dd2 = self.d["lens"].getLensDistance()
                if abs((Dd-Dd2)/Dd) > 1e-5 or abs((Dd-Dd2)/Dd2) > 1e-5:
                    raise LensInfoException("The angular diameter distance stored in the lens instance doesn't seem to match the specified 'zd'")

                Ds = self.cosm.getAngularDiameterDistance(self.d["zs"])
                Dds = self.cosm.getAngularDiameterDistance(self.d["zd"], self.d["zs"])
            else:
                raise LensInfoException("No distance info has been set; either Ds and Dds must be set, or zd and zs (and a cosmological model)")

            imgPlane = images.ImagePlane(lensPlane, Ds, Dds)
        else:
            zs = self.d["zs"]
            imgPlane = multiplane.MultiImagePlane(lensPlane, zs)

        self.d["imageplane"] = imgPlane
        return imgPlane

    def getIntegratedMass(self, renderer = "default", feedbackObject = "default"):
        """Based on the grid as set in the constructor, this returns the total
        estimated mass in the area, performing the necessary
        calculations if needed. With the `renderer` argument, you can specify a specific :ref:`renderer <renderers>`
        to speed up the calculation, and with `feedbackObject` you can specify a
        particular :ref:`feedback mechanism <feedback>`."""

        if "totalmass" in self.d:
            return self.d["totalmass"]

        Dd = self.d["lens"].getLensDistance()
        massPixelSum = sum(sum(self.getDensityPixels(renderer, feedbackObject)))

        bottomLeft = self.d["bottomleft"]
        topRight = self.d["topright"]
        totalWidth = abs(topRight[0] - bottomLeft[0])
        totalHeight = abs(topRight[1] - bottomLeft[1])

        self.d["totalmass"] = Dd**2 * (totalWidth/self.d["numx"])*(totalHeight/self.d["numy"]) * massPixelSum

        return self.d["totalmass"]

    def getDensityPoints(self, renderer = "default", feedbackObject = "default"):
        """Returns the densities on a 2D NumPy grid that corresponds to the
        area in the constructor, calculations if needed. With the `renderer` argument, 
        you can specify a specific :ref:`renderer <renderers>`
        to speed up the calculation, and with `feedbackObject` you can specify a
        particular :ref:`feedback mechanism <feedback>`."""

        if "densitypoints" in self.d:
            return self.d["densitypoints"]

        bottomLeft = self.d["bottomleft"]
        topRight = self.d["topright"]
        numX = self.d["numx"] + 1
        numY = self.d["numy"] + 1

        if not "lens" in self.d:
            raise LensInfoException("No lens has been set")

        lens = self.d["lens"]
        if type(self.d["lens"]) == list:
            raise LensInfoException("Can't calculate the density for a multiple lens plane lens")

        massMap = lens.getSurfaceMassDensityMap(bottomLeft, topRight, numX, numY, renderer = renderer, feedbackObject = feedbackObject, reduceToPixels = False)
        self.d["densitypoints"] = massMap

        return massMap

class PlotException(Exception):
    """An exception that will be thrown in case something goes wrong when plotting"""
    pass

# Based on https://stackoverflow.com/a/48695245/2828217
def plot3DInteractive(X, Y, Z, height=600, xlabel = "X", ylabel = "Y", zlabel = "Z", flipX = False, initialCamera = None, visJSoptions = None,
                      maxPoints = 128,
                      visJS="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.js",
                      visCSS="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.css",
                      canvas2svgJS="https://cdn.rawgit.com/gliffy/canvas2svg/master/canvas2svg.js",
                      fileSaverJS="https://cdn.rawgit.com/eligrey/FileSaver.js/b4a918669accb81f184c610d741a4a8e1306aa27/FileSaver.js"
                      ):
    """Helper function to creates an interative plot in a jupyter notebook, 
    using `vis.js <http://visjs.org/>`_.

    Arguments:
     - `X`: a NY x NX numpy array containing the X coordinates of the plot region.

     - `Y`: similar `X`, but with Y coordinates of the plot region.

     - `Z`: also a NY x NX array, containing the Z-values of the surface to be plotted.

     - `height`: the height of the plot area.

     - `xlabel`: the label for the x-axis.

     - `ylabel`: the label for the y-axis.

     - `zlabel`: the label for the z-axis.

     - `flipX`: a boolean that indicates if the x-axis should go from right to left.

     - `initialCamera`: when manipulating the view, you'll see information about 'horizontal',
       'vertical' and 'distance' information change. To make the initial view correspond to a
       particular setting, put these values in a dictionary with these keys.

     - `visJSoptions`: this is for advanced use only, you can pass extra options to a 
       `vis.Graph3d <http://visjs.org/docs/graph3d/#Configuration_Options>`_ instance this way.

     - `maxPoints`: an exception will be raised if the array dimensions exceed this value, this
       is to prevent accidentally plotting something with very high detail, which would become
       very slow.

     - `visJS`: URL to the vis.js JavaScript library.

     - `visCSS`: URL to the vis.js CSS file.

     - `canvas2svgJS`: URL to the canvas2svg.js file.

     - `fileSaverJS`: URL to FileSaver.js file
    """

    if X.shape[0] > maxPoints or X.shape[1] > maxPoints:
        raise PlotException("Number of points in the grid exceeds the maximum number of points; this is a warning to avoid the interactive view becoming too slow, override this (if you're certain!) by increasing the `maxPoints` parameter")

    if X.shape != Y.shape or X.shape != Z.shape:
        raise PlotException("The different arrays to be plotted do not have the same shape")

    options = {
        "width": "100%",
        "style": "surface",
        "showPerspective": True,
        "showGrid": True,
        "showShadow": False,
        "keepAspectRatio": True,
        "height": str(height) + "px",
        "xLabel": xlabel,
        "yLabel": ylabel,
        "zLabel": zlabel,
    }

    flipXCode = ""
    if flipX:
        flipXCode = """options["xValueLabel"] = function(x) 
        { 
            var s = "" + x;
            if (s[0] == "-")
                return s.substr(1);
            if (s == "0")
                return s;
            return "-" + s; 
        };""";

    if initialCamera:
        options["cameraPosition"] = initialCamera
        
    if visJSoptions:
        for x in visJSoptions:
            options[x] = copy.deepcopy(visJSoptions[x])
        
    flipFactor = -1.0 if flipX else 1.0
    dataJson = [ ]
    s = X.shape
    for y in range(s[0]):
        for x in range(s[1]):
            dataJson.append({"x": X[y,x]*flipFactor, "y": Y[y,x], "z": Z[y,x]})

    dataJson = json.dumps(dataJson)
    options = json.dumps(options)
    visCode = r"""
<html>
    <head>
        <link href=""" + '"' + visCSS + '"' + r""" type="text/css" rel="stylesheet" />
        <script src=""" + '"' + visJS + '"' + r""" type="text/javascript"></script>
        <script src=""" + '"' + canvas2svgJS + '"' + r""" type="text/javascript"></script>
        <script src=""" + '"' + fileSaverJS + '"' + r""" type="text/javascript"></script>
    </head>
    <body>
        <div id="pos" style="top:0px;left:0px;position:absolute;"></div>
        <div id="visualization"></div>
        <script type="text/javascript">
        var data = new vis.DataSet();
        data.add(""" + dataJson + r""");
        var options = """ + options + r""";

        """ + flipXCode + r"""

        var container = document.getElementById("visualization");
        graph3d = new vis.Graph3d(container, data, options);
        graph3d.on("cameraPositionChange", function(evt)
        {
            elem = document.getElementById("pos");
            s = "horizontal: " + evt.horizontal.toExponential(4) + "<br>vertical: " + evt.vertical.toExponential(4) + "<br>distance: " + evt.distance.toExponential(4);
            elem.innerHTML = s;
        });
        </script>
        <button onclick="exportSVG()" style="position:fixed;top:0px;right:0px;">Save to SVG</button>
        <script>
        function T(x)
        {
            var s = "" + x;
            while (s.length < 2)
                s = "0" + s;
            return s;
        }

        function exportSVG()
        {
            var cnvs = graph3d.frame.canvas;
            var fakeCtx = C2S(cnvs.width, cnvs.height);
            var realGetContext = cnvs.getContext;
            cnvs.getContext = function() { return fakeCtx; }

            graph3d.redraw();
            var svg = fakeCtx.getSerializedSvg();

            cnvs.getContext = realGetContext;
            graph3d.redraw();

            var b = new Blob([svg], { type: "image/svg+xml;charset=utf-8" });
            var d = new Date();
            var fileName = "Capture-" + d.getFullYear() + "-" + T(d.getMonth()+1) + "-" + T(d.getDate()) + "_" + T(d.getHours()) + "-" + T(d.getMinutes()) + "-" + T(d.getSeconds()) + ".svg";
            saveAs(b, fileName);
        }
        </script>
    </body>
</html>"""
    #print(visCode)
    htmlCode = "<iframe srcdoc='"+visCode+"' width='100%' height='" + str(height) + "px' style='border:0;' scrolling='no'> </iframe>"

    from IPython.core.display import display, HTML
    display(HTML(htmlCode))

def plotDensityInteractive(lensOrLensInfo, numX=75, numY=75, height=600, xlabel="X", ylabel="Y", zlabel="Z", 
        angularUnit="default", densityUnit=1.0, renderer="default", feedbackObject="default", flipX = False,
        initialCamera = None, visJSoptions = None, maxPoints = 128, 
        visJS="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.js",
        visCSS="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.css",
        canvas2svgJS="https://cdn.rawgit.com/gliffy/canvas2svg/master/canvas2svg.js",
        fileSaverJS="https://cdn.rawgit.com/eligrey/FileSaver.js/b4a918669accb81f184c610d741a4a8e1306aa27/FileSaver.js"
        ):
    """Creates an interactive 3D plot of the mass density specified in `lensOrLensInfo`.

    Arguments:
     - `lensOrLensInfo`: this can either be a :class:`GravitationalLens <grale.lenses.GravitationalLens>`
       instance or an instance of :class:`LensInfo` or even :class:`DensInfo`.
       In case it's only a gravitational lens, an estimate of the relevant area will be used.

     - `numX`: the number of points in the x-direction used to create the surface or 2D
       plot.

     - `numY`: the number of points in the y-direction used to create the surface or 2D
       plot.

     - `height`: the height of the plot area.

     - `xlabel`: the label for the x-axis.

     - `ylabel`: the label for the y-axis.

     - `zlabel`: the label for the z-axis.

     - `angularUnit`: the angular unit that should be used in the plot. The 
       :ref:`pre-defined constants <constants>` can be useful here.

     - `densityUnit`: by default, the density will be in kg/m^2, but another unit can be specified
       here.

     - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
       to speed up the calculation. 

     - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

     - `flipX`: a boolean that indicates if the x-axis should go from right to left.

     - `initialCamera`: when manipulating the view, you'll see information about 'horizontal',
       'vertical' and 'distance' information change. To make the initial view correspond to a
       particular setting, put these values in a dictionary with these keys.

     - `visJSoptions`: this is for advanced use only, you can pass extra options to a 
       `vis.Graph3d <http://visjs.org/docs/graph3d/#Configuration_Options>`_ instance this way.

     - `maxPoints`: an exception will be raised if the array dimensions exceed this value, this
       is to prevent accidentally plotting something with very high detail, which would become
       very slow.

     - `visJS`: URL to the vis.js JavaScript library.

     - `visCSS`: URL to the vis.js CSS file.

     - `canvas2svgJS`: URL to the canvas2svg.js file.

     - `fileSaverJS`: URL to FileSaver.js file
    """
    angularUnit = _getAngularUnit(angularUnit)

    lensInfo = _toLensInfo(lensOrLensInfo)

    # TODO: can we replace this using a gridfunction?
    plotArray = privutilcython.resample2DArray(lensInfo.getDensityPoints(renderer, feedbackObject), numY, numX)
    bottomLeft = lensInfo.getBottomLeft()
    topRight = lensInfo.getTopRight()

    X, Y = np.meshgrid(np.linspace(bottomLeft[0], topRight[0], numX),
                       np.linspace(bottomLeft[1], topRight[1], numY))

    X /= angularUnit
    Y /= angularUnit
    Z = plotArray/densityUnit
    
    plot3DInteractive(X, Y, Z, flipX=flipX, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, initialCamera=initialCamera, 
                      visJSoptions=visJSoptions, maxPoints=maxPoints, height=height, visJS=visJS, visCSS=visCSS,
                      canvas2svgJS=canvas2svgJS, fileSaverJS=fileSaverJS)
    return lensInfo

def plotDensityContours(lensOrLensInfo, renderer = "default", feedbackObject = "default", angularUnit = "default", densityUnit = 1.0,
                axes = None, axImgCallback = None, levels = None, **kwargs):
    """Creates a contour plot of the mass density specified in `lensOrLensInfo`.

In case you just want the calculations to be performed (for example because you need
the calculated density pixels) but you don't want an actual plot, you can set the
`axes` parameter to `False`.

Arguments:
 - `lensOrLensInfo`: this can either be a :class:`GravitationalLens <grale.lenses.GravitationalLens>`
   instance or an instance of :class:`LensInfo` or even :class:`DensInfo`.
   In case it's only a gravitational lens, an estimate of the relevant area will be used.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `densityUnit`: by default, the density will be in kg/m^2, but another unit can be specified
   here.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. The value `False` has a special meaning: in that case,
   the calculations will be performed as usual, but an actual plot will not be created.

 - `axImgCallback`: if specified, this callback function will be called with the object
   returned by ``contour`` as argument.

 - `kwargs`: these parameters will be passed on to the `contour <https://matplotlib.org/api/_as_gen/matplotlib.pyplot.contour.html>`_
   function in matplotlib.
"""
    angularUnit = _getAngularUnit(angularUnit)

    lensInfo = _toLensInfo(lensOrLensInfo)

    if axes is not False:

        X = lensInfo.getXPointCoordinates()/angularUnit
        Y = lensInfo.getYPointCoordinates()/angularUnit
        Z = lensInfo.getDensityPoints(renderer, feedbackObject)/densityUnit

        bottomLeft = lensInfo.getBottomLeft()
        topRight = lensInfo.getTopRight()

        if axes is None:
            import matplotlib.pyplot as plt
            axes = plt.gca()

        # Note: need to swap Y labeling here because of the way the pixels are ordered in this function
        if levels is None:
            axImg = axes.contour(X, Y, Z, **kwargs)
        else:
            axImg = axes.contour(X, Y, Z, levels, **kwargs)

        axes.set_xlim([bottomLeft[0]/angularUnit, topRight[0]/angularUnit])
        axes.set_ylim([bottomLeft[1]/angularUnit, topRight[1]/angularUnit])
        axes.set_aspect("equal")

        if axImgCallback: axImgCallback(axImg)

    return lensInfo

def plotDensity(lensOrLensInfo, renderer = "default", feedbackObject = "default", angularUnit = "default", densityUnit = 1.0,
                axes = None, axImgCallback = None, **kwargs):
    """Creates a 2D plot of the situation specified in `lensOrLensInfo`.

In case you just want the calculations to be performed (for example because you need
the calculated density pixels) but you don't want an actual plot, you can set the
`axes` parameter to `False`.

Arguments:
 - `lensOrLensInfo`: this can either be a :class:`GravitationalLens <grale.lenses.GravitationalLens>`
   instance or an instance of :class:`LensInfo` or even :class:`DensInfo`.
   In case it's only a gravitational lens, an estimate of the relevant area will be used.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `densityUnit`: by default, the density will be in kg/m^2, but another unit can be specified
   here.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. The value `False` has a special meaning: in that case,
   the calculations will be performed as usual, but an actual plot will not be created.

 - `axImgCallback`: if specified, this callback function will be called with the object
   returned by ``imshow`` as argument.

 - `kwargs`: these parameters will be passed on to the `imshow <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.imshow.html>`_
   function in matplotlib.
"""
    angularUnit = _getAngularUnit(angularUnit)

    lensInfo = _toLensInfo(lensOrLensInfo)

    if axes is not False:

        pixels = lensInfo.getDensityPixels(renderer, feedbackObject)
        bottomLeft = lensInfo.getBottomLeft()
        topRight = lensInfo.getTopRight()

        if axes is None:
            import matplotlib.pyplot as plt
            axes = plt.gca()

        # Note: need to swap Y labeling here because of the way the pixels are ordered in this function
        axImg = axes.imshow(pixels/densityUnit, extent = np.array([ bottomLeft[0], topRight[0], topRight[1], bottomLeft[1]])/angularUnit, **kwargs)
        axes.set_xlim([bottomLeft[0]/angularUnit, topRight[0]/angularUnit])
        axes.set_ylim([bottomLeft[1]/angularUnit, topRight[1]/angularUnit])   
        axes.set_aspect("equal")
        if axImgCallback: axImgCallback(axImg)

    return lensInfo

def plotImagePlane(lensOrLensInfo, sources = [], renderer = "default", feedbackObject = "default", 
                   angularUnit = "default", subSamples = 9, sourceRgb = (0, 1, 0), imageRgb = (1, 1, 1),
                   bgRgb = (0, 0, 0), caustColor = "blue", caustKw = {}, critColor = "red", critKw = {},
                   plotCaustics = True, plotCriticalLines = True, plotSources = True, plotImages = True,
                   evenError = True, axes = None, axImgCallback = None, 
                   processRenderPixels = None, **kwargs):
    """Create a matplotlib-based plot of image plane and/or source plane for certain
lens parameters in `lensOrLensInfo`. You can also use this function to create the necessary
calculated mappings but not the plot, by setting `axes` to `False`. This can be useful
to create a :py:class:`LensPlane<grale.images.LensPlane>` and :py:class:`ImagePlane<grale.images.ImagePlane>` instance
using a specific renderer to speed up the calculation.

Arguments:
 - `lensOrLensInfo`: this can either be a :class:`GravitationalLens <grale.lenses.GravitationalLens>`
   instance or an instance of :class:`LensInfo`.
   In case it's only a gravitational lens, an estimate of the corners will be used, and
   both `Ds` and `Dds` will be set to 1.0 (only the value of Dds/Ds matters, and will be
   one in this case).

 - `sources`: one or more sources for which the lens effect will be calculated.
   A source can either be a :ref:`source shape <sourceshapes>`, or a dictionary with
   the following entries:

       - ``shape``: the actual :ref:`source shape <sourceshapes>`
       - ``z``, or ``Ds`` and  ``Dds``: redshift or angular diameter distances for this
         source.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `subSamples`: each pixel is sub-sampled sqrt(subSamples) times in x- and y- direction, to be able to 
   roughly approximate the integration that’s needed over the surface area of a pixel. This parameter
   is passed to :py:meth:`ImagePlane.renderImages<grale.images.ImagePlane.renderImages>`
   and :py:meth:`ImagePlane.renderSources<grale.images.ImagePlane.renderSources>`.

 - `sourceRgb`: the RGB color, or a list of RGB colors, (each component is a number between 0 and 1) 
   that should be given to a pixel that lies within a source.

 - `imageRgb`: the RGB color, or a list of RGB colors, (each component is a number between 0 and 1) 
   that should be given to a pixel that lies within an image.
 
 - `caustColor`: color, or list of colors, for the caustics.

 - `caustKw`: dictionary describing keyword arguments for the matplotlib's ``plot`` function
   for drawing the caustics.

 - `critColor`: color, or list of colors, for the critical lines.

 - `critKw`: dictionary describing keyword arguments for the matplotlib's ``plot`` function
   for drawing the critical lines.

 - `bgRgb`: the RGB color for the background.

 - `plotCaustics`: boolean value, or list of values, that indicates if the caustics 
   for a source should be drawn on the plot.

 - `plotCriticalLines`: boolean value, or list of values, that indicates if the critical lines 
   for a source should be drawn on the plot.

 - `plotSources`: boolean value (or list) that indicates if the sources should be drawn on the plot.

 - `plotImages`: boolean value (or list) that indicates if the images should be drawn on the plot.

 - `evenError`: by default, the function will raise an exception if the number of pixels
   in the x or y-direction is even, because for simple lenses odd pixel numbers work better.
   In case you don't care for this, set the flag to `False` and this will stop the calculation.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. The value `False` has a special meaning: in that case,
   the calculations will be performed as usual, but an actual plot will not be created.

 - `axImgCallback`: if specified, this callback function will be called with the object
   returned by ``imshow`` as argument.

 - `processRenderPixels`: a function that is called before the final ``imshow`` call; it
   is the array that's returned by this function (if specified) that is plotted. This can
   be used to combine the images of several sources with a galaxy image for example.

 - `kwargs`: these parameters will be passed on to the `imshow <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.imshow.html>`_
   function in matplotlib.
"""
    angularUnit = _getAngularUnit(angularUnit)

    lensInfo = _toLensInfo(lensOrLensInfo)

    bottomLeft = lensInfo.getBottomLeft()
    topRight = lensInfo.getTopRight()
    numX = lensInfo.getNumXPixels()
    numY = lensInfo.getNumYPixels()

    if axes is None:
        import matplotlib.pyplot as plt
        axes = plt.gca()

    # TODO: sourceScale and imageScale
    
    f = lambda x : np.array(x, dtype=np.double)
    bgRgb = f(bgRgb)
    
    if len(bgRgb) == 1:
        bgPlane = np.multiply(np.ones((numY, numX)), bgRgb)
    else:
        bgPlane = np.multiply(np.ones((numY, numX, 1)), bgRgb)

    def setSrcDistAndGetImagePlane(srcEntry):
        if "z" in srcEntry:
            lensInfo.setSourceRedshift(srcEntry["z"])
        elif "Dds" in srcEntry and "Ds" in srcEntry:
            lensInfo.setSourceDistances(srcEntry["Ds"], srcEntry["Dds"])
        else:
            raise PlotException("Expecting either 'z' or 'Ds' and 'Dds' to be present in the source information")
        
        return lensInfo.getImagePlane(renderer, feedbackObject)

    if not sources:
        # Make sure imgPlane exists
        imgPlane = lensInfo.getImagePlane(renderer, feedbackObject)
    else:
        def createRgbPlane(plane, rgb):
            if len(rgb) != 1:
                plane = plane.reshape(plane.shape + (1,))

            return np.multiply(plane, rgb)

        if type(sources) != list:
            sources = [ sources ]

        if type(sources[0]) == dict:

            if len(bgRgb) == 1:
                rgbSplane = np.zeros((numY, numX))
                rgbIplane = np.zeros((numY, numX))
            else:
                rgbSplane = np.zeros((numY, numX, len(bgRgb)))
                rgbIplane = np.zeros((numY, numX, len(bgRgb)))

            sourceRgbList = sourceRgb if type(sourceRgb) == list else [ sourceRgb for s in sources ]
            imgRgbList = imageRgb  if type(imageRgb) == list else [ imageRgb for s in sources ]

            for srcEntry, sourceRgb, imageRgb in zip(sources, sourceRgbList, imgRgbList):
                splane, iplane = [ None ], [ None ]
                imgPlane = setSrcDistAndGetImagePlane(srcEntry)
                
                for flag, renderFunction, destPlaneList in [
                        (plotSources, imgPlane.renderSources, splane),
                        (plotImages, imgPlane.renderImages, iplane)
                        ]:
                    if flag:
                        plane = renderFunction([srcEntry["shape"]], subSamples = subSamples)
                        if destPlaneList[0] is None:
                            destPlaneList[0] = plane
                        else:
                            destPlaneList[0] += plane

                splane = splane[0]
                iplane = iplane[0]
                                
                if splane is not None:
                    rgbSplane += createRgbPlane(splane, sourceRgb)
                if iplane is not None:
                    rgbIplane += createRgbPlane(iplane, imageRgb)

        else:
            imgPlane = lensInfo.getImagePlane(renderer, feedbackObject)

            if plotSources:
                splane = imgPlane.renderSources(sources, subSamples = subSamples)
                rgbSplane = createRgbPlane(splane, sourceRgb)

            if plotImages:
                iplane = imgPlane.renderImages(sources, subSamples = subSamples)
                rgbIplane = createRgbPlane(iplane, imageRgb)

    if axes is not False:
        # This allows information from multiple planes to be accumulated
        plane = processRenderPixels(plane) if processRenderPixels else np.clip(bgPlane+rgbSplane+rgbIplane, 0, 1)
        # Note: need to swap Y labeling here because of the way the pixels are ordered in this function
        if plane is not None:
            axImg = axes.imshow(plane, extent = np.array([ bottomLeft[0], topRight[0], topRight[1], bottomLeft[1]])/angularUnit, **kwargs)
            if axImgCallback: axImgCallback(axImg)
        
    def plotLines(lines, angularUnit, **kwargs):
        for part in lines:
            x, y = [], []
            for a in part:
                x.append(a[0]/angularUnit)
                y.append(a[1]/angularUnit)
#                print(x[-1],y[-1])

            if axes is not False:
                axes.plot(x, y, **kwargs)
#            print()

    for flagInfo, isCaust, colorInfo, keywords in [
            (plotCaustics, True, caustColor, caustKw),
            (plotCriticalLines, False, critColor, critKw)
            ]:
        if flagInfo:
            if sources and type(sources[0]) == dict:
                if type(flagInfo) != list:
                    flagInfo = [ flagInfo for s in sources ]
                if type(colorInfo) != list:
                    colorInfo = [ colorInfo for c in colorInfo ]

                for srcEntry, plotFlag, col in zip(sources, flagInfo, colorInfo):
                    if plotFlag:
                        imgPlane = setSrcDistAndGetImagePlane(srcEntry)
                        lines = imgPlane.getCaustics() if isCaust else imgPlane.getCriticalLines()
                        plotLines(lines, angularUnit, color=col, **keywords)

            else:
                lines = imgPlane.getCaustics() if isCaust else imgPlane.getCriticalLines()
                plotLines(lines, angularUnit, color=colorInfo, **keywords)


    if axes is not False:
        axes.set_xlim([bottomLeft[0]/angularUnit, topRight[0]/angularUnit])
        axes.set_ylim([bottomLeft[1]/angularUnit, topRight[1]/angularUnit])   
        axes.set_aspect("equal")

    return lensInfo

def _getGnuplotVersion(gnuplotExe):
    p = subprocess.Popen([gnuplotExe , "--version"], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    v = S(p.communicate()[0])
    version = int(v.splitlines()[0].split()[1].split(".")[0])
    return version

def plotImagePlaneGnuplot(lensOrLensInfo, fileNameBase = None, sources = [], 
                          plotCaustics = True, plotCriticalLines = True, plotSources = True, plotImages = True,
                          renderer = "default", feedbackObject = "default", angularUnit = "default", subSamples = 9,
                          evenError = True, axes = None, color = True, xlabel = "X", ylabel = "Y", grid = False,
                          flipX = False, gnuplotExe = "gnuplot", convertExe = "convert"):

    """Creates a gnuplot-based plot of the image plane/source plane that corresponds to the
situation specified in `lensOrLensInfo`, saving various files with names starting with 
`fileNameBase`. If you only want to generate these files, but not visualize them in a plot 
immediately, you can set `axes` to `False`.

Arguments:
 - `lensOrLensInfo`: this can either be a :class:`GravitationalLens <grale.lenses.GravitationalLens>`
   instance or an instance of :class:`LensInfo`.
   In case it's only a gravitational lens, an estimate of the corners will be used, and
   both `Ds` and `Dds` will be set to 1.0 (only the value of Dds/Ds matters, and will be
   one in this case).

 - `fileNameBase`: the generated gnuplot file will have this name, ending with ``.gnuplot``.
   Similarly, the generated eps file will start with this and end with ``.eps``. If set to
   `None`, no output files will be generated, but the result may still be shown.

 - `sources`: a list of :ref:`source shapes <sourceshapes>` that should be used to calculate
   the images from.

 - `plotCaustics`: boolean value that indicates if the caustics should be drawn on the plot.

 - `plotCriticalLines`: boolean value that indicates if the critical lines should be drawn on
   the plot.

 - `plotSources`: boolean value that indicates if the sources should be drawn on the plot.

 - `plotImages`: boolean value that indicates if the images should be drawn on the plot.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `subSamples`: each pixel is sub-sampled sqrt(subSamples) times in x- and y- direction, to be able to 
   roughly approximate the integration that’s needed over the surface area of a pixel. This parameter
   is passed to :py:meth:`ImagePlane.renderImages<grale.images.ImagePlane.renderImages>`
   and :py:meth:`ImagePlane.renderSources<grale.images.ImagePlane.renderSources>`.

 - `evenError`: by default, the function will raise an exception if the number of pixels
   in the x or y-direction is even, because for simple lenses odd pixel numbers work better.
   In case you don't care for this, set the flag to `False` and this will stop the calculation.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. The value `False` has a special meaning: in that case,
   the calculations will be performed as usual, but an actual plot will not be created.

 - `color`: a boolean that indicates if the plot should use colors or grayscale only.

 - `xlabel`: the label for the x-axis.

 - `ylabel`: the label for the y-axis.

 - `grid`: a boolean that indicates if a grid should be shown as well.

 - `flipX`: a boolean that indicates if the x-axis should go from right to left.

 - `gnuplotExe`: the ``gnuplot`` command.

 - `convertExe`: the ``convert`` command. This program is used to show the generated plot
   using matplotlib: the eps generated by gnuplot is converted to a PNG file which is then
   shown using matplotlib.

"""
    angularUnit = _getAngularUnit(angularUnit)

    lensInfo = _toLensInfo(lensOrLensInfo)
    gpVersion = _getGnuplotVersion(gnuplotExe)

    bottomLeft = lensInfo.getBottomLeft()
    topRight = lensInfo.getTopRight()
    imgPlane = lensInfo.getImagePlane(renderer, feedbackObject)
    renderInfo = imgPlane.getRenderInfo()

    pointsize = 70.0/renderInfo["xpoints"]
    if pointsize < 0.5:
        pointsize = 0.5

    fontsize = 14
    pointtype = 31

    gnuplotData = ""
    if color:
        gnuplotData += "set terminal postscript eps color solid enhanced {}\n".format(fontsize)
    else:
        gnuplotData += "set terminal postscript eps monochrome enhanced {}\n".format(fontsize)

    gnuplotData += """
set size 1.0,1.3
set style data points
set pointsize {}
set xlabel '{}'
set ylabel '{}'
""".format(pointsize, xlabel, ylabel)

    if grid:
        gnuplotData += "set grid\n"
    
    if gpVersion >= 5:
        gnuplotData += "set colors classic\n"
        lw = 1.5
    else:
        lw = 3

    first = True
    count = 0
    for i in [ plotCaustics, plotCriticalLines, plotSources, plotImages ] : 
        if i:
            count += 1

    x1, x0 = topRight[0]/angularUnit, bottomLeft[0]/angularUnit
    if flipX:
        x0, x1 = x1, x0

    y1, y0 = topRight[1]/angularUnit, bottomLeft[1]/angularUnit

    gnuplotData += "plot [{}:{}][{}:{}] ".format(x0,x1,y0,y1)

    if count == 0:
        gnuplotData += "{} notitle\n".format(y0 - (y1-y0))
        return

    if plotCaustics:
        first = False
        if color:
            gnuplotData += "'-' t 'Caustics' with linespoints lw {} lt 3 pt 1 ps 0.05".format(lw)
        else:
            gnuplotData += "'-' t 'Caustics' with linespoints lw {} lt 1 pt 1 ps 0.05".format(lw)

    if plotCriticalLines:
        if not first: gnuplotData += ","
        first = False
        
        if color:
            gnuplotData += "'-' t 'Critical lines' with linespoints lw {} lt 1 pt 1 ps 0.05".format(lw)
        else:
            gnuplotData += "'-' t 'Critical lines' with lines lw {} lt 3".format(lw)

    if plotSources:
        if not first: gnuplotData += ","
        first = False

        if color:
            gnuplotData += "'-' t 'Sources' pt {} lt 9".format(pointtype)
        else:
            gnuplotData += "'-' t 'Sources' pt {} lt 1".format(pointtype)

    if plotImages:
        if not first: gnuplotData += ","
        first = False

        if color:
            gnuplotData += "'-' t 'Images' pt {} lt 0".format(pointtype)
        else:
            gnuplotData += "'-' t 'Images' pt {} lt 1".format(pointtype)

    gnuplotData += "\n"

    def plotCritCaust(lines):
        gnuplotData = ""
        first = True
        for line in lines:
            if not first: gnuplotData += "\n"
            first = False
            for p in line:
                gnuplotData += "{}\t{}\n".format(float(p[0])/angularUnit, float(p[1])/angularUnit)

        if len(lines) == 0:
            gnuplotData += "{}\t{}\n".format(float(topRight[0])/angularUnit*2, float(topRight[1])/angularUnit*2)

        gnuplotData += "e\n"
        return gnuplotData

    if plotCaustics:
        gnuplotData += plotCritCaust(imgPlane.getCaustics())

    if plotCriticalLines:
        gnuplotData += plotCritCaust(imgPlane.getCriticalLines())

    splane = imgPlane.renderSources(sources, subSamples = subSamples)
    iplane = imgPlane.renderImages(sources, subSamples = subSamples)

    def plotSourcesImages(plane):

        gnuplotData = ""

        numy = plane.shape[0]
        numx = plane.shape[1]

        pixw = (topRight[0]-bottomLeft[0])/numx
        pixh = (topRight[1]-bottomLeft[1])/numy

        yis, xis = np.where(plane > 0)
        for i in range(len(yis)):
            xi = xis[i]
            yi = yis[i]
            
            x = (bottomLeft[0] + (0.5+xi)*pixw)/angularUnit
            y = (bottomLeft[1] + (0.5+yi)*pixh)/angularUnit
            v = plane[yi,xi]

            gnuplotData += "{}\t{}\t{}\n".format(x, y, v)

        if len(yis) == 0:
            gnuplotData += "{}\t{}\t0\n".format(topRight[0]/angularUnit*2, topRight[1]/angularUnit*2)

        gnuplotData += "e\n"
        return gnuplotData

    if plotSources:
        gnuplotData += plotSourcesImages(splane)

    if plotImages:
        gnuplotData += plotSourcesImages(iplane)

    return (lensInfo, ) + _finalizeGnuplotPlot(fileNameBase, gnuplotData, axes, gnuplotExe, convertExe)

def _finalizeGnuplotPlot(fileNameBase, gnuplotData, axes, gnuplotExe, convertExe, rotateAngle = 0):

    filesToDelete = [ ]
    noOutput = False

    if fileNameBase is None:
        noOutput = True
        f = tempfile.NamedTemporaryFile("w+t", delete = False)
        fileNameBase = f.name
        f.close()

    gnuplotFileName = fileNameBase + ".gnuplot"
    epsFileName = fileNameBase + ".eps"
    pngFileName = fileNameBase + ".png"

    if noOutput:
        filesToDelete.append(fileNameBase)
        filesToDelete.append(gnuplotFileName)
        filesToDelete.append(epsFileName)
        # filesToDelete.append(pngFileName) # This will be deleted anyway

    with open(gnuplotFileName, "wt") as f:
        f.write(gnuplotData)

    with open(epsFileName, "wt") as f:
        subprocess.check_call([gnuplotExe, gnuplotFileName], stdout = f)

    if axes is not False:
        try:
            subprocess.check_call([convertExe, "-flatten", "-density", "200x200", epsFileName, "-rotate", str(rotateAngle), "-define", "png:color-type=2", pngFileName])
            subprocess.check_call([convertExe, "-trim", pngFileName,  "-define", "png:color-type=2",pngFileName])

            import matplotlib.pyplot as plt
            if axes is None:
                axes = plt.gca()

            axes.imshow(plt.imread(pngFileName))
            axes.axis('off')

            os.unlink(pngFileName)
        except Exception as e:
            print("Warning: unable to visualize results directly: {}".format(e))

    for n in filesToDelete:
        if os.path.exists(n):
            try:
                os.unlink(n)
            except Exception as e:
                print("Warning, unable to remove file '{}': {}".format(n, e))

    if noOutput:
        return (None, None)

    return (gnuplotFileName, epsFileName) 

def plotDensityGnuplot(lensOrLensInfo, fileNameBase = None, numX = 75, numY = 75, renderer = "default", 
                       feedbackObject = "default", angularUnit = "default", densityUnit = 1.0,
                       axes = None, color = True, xlabel = "X", ylabel = "Y", zlabel = "Z",
                       plot3D = True, rotationX = 45, rotationZ = 45, showAxes = True,
                       flipX = False, gnuplotExe = "gnuplot", convertExe = "convert"):
    """Creates a gnuplot-based plot of the mass density of the lens, either as a 3D surface
plot or as a 2D plot. Various files with names starting with `fileNameBase` are
generated. If you only want to generate these files, but not visualize them in a plot 
immediately, you can set `axes` to `False`.

Arguments:
 - `lensOrLensInfo`: this can either be a :class:`GravitationalLens <grale.lenses.GravitationalLens>`
   instance or an instance of :class:`LensInfo` or even :class:`DensInfo`.
   In case it's only a gravitational lens, an estimate of the relevant area will be used.

 - `fileNameBase`: the generated gnuplot file will have this name, ending with ``.gnuplot``.
   Similarly, the generated eps file will start with this and end with ``.eps``. If set to
   `None`, no output files will be generated, but the result may still be shown.

 - `numX`: the number of points in the x-direction used to create the surface or 2D
   plot.

 - `numY`: the number of points in the y-direction used to create the surface or 2D
   plot.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `densityUnit`: by default, the density will be in kg/m^2, but another unit can be specified
   here.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. The value `False` has a special meaning: in that case,
   the calculations will be performed as usual, but an actual plot will not be created.

 - `color`: a boolean that indicates if the plot should use colors or grayscale only.

 - `xlabel`: the label for the x-axis.

 - `ylabel`: the label for the y-axis.

 - `zlabel`: the label for the z-axis.

 - `plot3D`: flag indicating if a 3D surface plot should be generated, or a 2D plot.

 - `rotationX`: if a 3D surface plot is created, this specifies the number of degrees
   the plot is rotated around the x-axis before showing.

 - `rotationZ`: if a 3D surface plot is created, this specifies the number of degrees
   the plot is rotated around the z-axis before showing.

 - `showAxes`: flag that indicates if the axes should be shown on the plot.

 - `flipX`: a boolean that indicates if the x-axis should go from right to left.

 - `gnuplotExe`: the ``gnuplot`` command.

 - `convertExe`: the ``convert`` command. This program is used to show the generated plot
   using matplotlib: the eps generated by gnuplot is converted to a PNG file which is then
   shown using matplotlib.
"""
    angularUnit = _getAngularUnit(angularUnit)

    lensInfo = _toLensInfo(lensOrLensInfo)
    gpVersion = _getGnuplotVersion(gnuplotExe)
    
    # TODO: can we replace this using a gridfunction?
    plotArray = privutilcython.resample2DArray(lensInfo.getDensityPoints(renderer, feedbackObject), numY, numX)
    bottomLeft = lensInfo.getBottomLeft()
    topRight = lensInfo.getTopRight()

    pointsize = 75.0/numX
    gnuplotData = ""
    if not color:
        gnuplotData += "set terminal postscript enhanced monochrome solid\n"
    else:
        gnuplotData += "set terminal postscript enhanced color solid\n"

    gnuplotData += "set size 1.0,1.3\n"

    if gpVersion >= 5:
        gnuplotData += "set colors classic\n"

    if not plot3D:
        gnuplotData += "set pointsize {}\nset view map\n".format(pointsize)
    else:
        gnuplotData += """set hidden3d
set view {}, {}
set ticslevel 0
""".format(rotationX, rotationZ)

    if not showAxes:
        gnuplotData += """unset xtics
unset ytics
unset ztics
"""

    gnuplotData += "set noborder\n"
    if not color:
        gnuplotData += "set palette defined ( 0 'white',1 'black')\n"

    if xlabel:
        gnuplotData += "set xlabel '{}'\n".format(xlabel)
    if ylabel:
        gnuplotData += "set ylabel '{}'\n".format(ylabel)
    if zlabel:
        gnuplotData += "set zlabel '{}'\n".format(zlabel)

    if not plot3D:
        gnuplotData += "splot '-' notitle with points pt 47 palette\n"
    else:
        colSuffix = "palette" if color else ""
        xRange = "{}:{}".format(bottomLeft[0]/angularUnit, topRight[0]/angularUnit) if not flipX else "{}:{}".format(topRight[0]/angularUnit, bottomLeft[0]/angularUnit)
        yRange = "{}:{}".format(bottomLeft[1]/angularUnit, topRight[1]/angularUnit)
        zRange = "" # Make this configurable?

        gnuplotData += "splot [{}][{}][{}] '-' notitle with lines {}\n".format(xRange, yRange, zRange, colSuffix)

    for yi in range(numY):
        y = bottomLeft[1] + (topRight[1]-bottomLeft[1])*float(yi)/(numY-1)
        y /= angularUnit
        for xi in range(numX):
            x = bottomLeft[0] + (topRight[0]-bottomLeft[0])*float(xi)/(numX-1)
            x /= angularUnit

            z = plotArray[yi, xi]/densityUnit

            gnuplotData += "{} {} {}\n".format(x, y, z)

        gnuplotData += "\n"

    return (lensInfo,) + _finalizeGnuplotPlot(fileNameBase, gnuplotData, axes, gnuplotExe, convertExe, 90)

class Animation(object):
    """A helper class to create an animation. This class generates an mp4 or webm file
as soon as an instance is created. To visualize the result immediately in an notebook,
you can use the :py:class:`NotebookAnimation<grale.plotutil.NotebookAnimation>` class.

To generate a useful plot, you'll need to override at least the `onDraw` member function.
Additional initialization can be done in the `onInit` member function.

While the class can be used to generated animations that involve gravitational lensing,
this is certainly not required. Example:

.. code-block:: python

	import grale.plotutil as plotutil
	import numpy as np

	class MyAnim(grale.plotutil.Animation):
		def __init__(self):
			super(MyAnim, self).__init__("test.mp4", [0], [10], 100, 25)
		
		def onInit(self, axes):
			pass
		
		def onDraw(self, pos, axes):
			x = np.linspace(0, 6, 100)
			y = np.cos(x+pos)
			axes.clear()
			axes.plot(x, y)

	# Constructing an instance of this class will immediately start rendering the frames.
	anim = MyAnim()
"""
    def __init__(self, fileName, startPos, endPos, steps, fps):
        """Initializes the object and immediately starts rendering the frames into a movie.

Arguments:
 - `fileName`: the name of the file that will contain the animation. By default an mp4
   movie will be created, but if the filename ends with ``.webm``, a webm format will
   be used instead.

 - `startPos`, `endpos` and `steps`: for each frame, a position that goes from `startpos`
   to `endpos` in the number of iterations specified by `steps`, will be passed to the
   `onDraw` member function. Both `startpos` and `endpos` should be arrays, tuples or
   vectors, with the same dimension.

 - `fps`: the frames generated by the animation (will be `steps` frames) will be played
   back at this number of frames per second.
"""
        self.startPos = startPos
        self.endPos = endPos
        self.steps = steps
        
        import matplotlib.pyplot as plt
        self.fig = plt.gcf()
        self.axes = plt.gca()

        import matplotlib.animation as animation
        anim = animation.FuncAnimation(self.fig, self._animationStep, init_func=self._onInitInternal, frames=self.steps, interval = 1, blit = False)
        plt.close(anim._fig)

        extraArgs = [ '-vcodec', 'libx264' ]
        self.mimeType = 'video/mp4'
        if fileName.lower().endswith(".webm"):
            extraArgs = ['-vcodec', 'vp8']
            self.mimeType = 'video/webm'

        with open(fileName, "wb") as f:
            anim.save(f.name, fps = fps, extra_args=extraArgs)

    def _onInitInternal(self):
        self.onInit(self.axes)
        return self.fig,

    def _animationStep(self, i):
        positionObject = [ 0 ] * len(self.startPos)
        for j in range(len(positionObject)):
            x0 = self.startPos[j]
            x1 = self.endPos[j]
            positionObject[j] = float(i)/(self.steps-1) * (x1-x0) + x0

        self.onFrame(i, self.steps)
        self.onDraw(positionObject, self.axes)
        return self.fig,

    def onInit(self, axes):
        """In case some extra initialization needs to be done, it can happen in this
member function. This is called before the first `onDraw` call. The `axes` object is
the matplotlib `axes <https://matplotlib.org/api/axes_api.html>`_ object that needs
to be drawn on in the `onDraw` function."""
        pass

    def onDraw(self, pos, axes):
        """For each frame, a position between `startPos` and `endPos` (from the constuctor)
is determined and this value is passed as `pos`. You can then draw something on the 
`axes` object that corresponds to this position."""
        pass

    def onFrame(self, frameNum, totalFrames):
        """This method is called to provide you with feedback. By default it writes the
progress to stdout; override to do something different."""
        print("{}/{}".format(frameNum+1,totalFrames))

class NotebookAnimation(Animation):
    """This is derived from the :py:class:`Animation <grale.plotutil.Animation>` class
and not only creates the movie but also embeds it in a notebook for immediate visualization.
The constructor takes the same arguments as the one from :py:class:`Animation <grale.plotutil.Animation>`.

Note that only a reference to the movie is embedded, so the movie will not be saved inside
the notebook itself.
"""
    def __init__(self, fileName, *args, **kwargs):
        from IPython.display import HTML, display

        self.feedback = None
        super(NotebookAnimation, self).__init__(fileName, *args, **kwargs)
        self.feedback = None # This should clear the progress bar and status field

        display(HTML("<video controls src='{}' type='{}'></video>".format(fileName, self.mimeType)))

    def onFrame(self, frameNum, totalFrames):
        if frameNum == 0:
            self.feedback = feedback.NotebookFeedback(1, totalFrames)

        self.feedback.onProgress(frameNum+1)
        self.feedback.onStatus("{}/{}".format(frameNum+1,totalFrames))

def estimatePlotScale(lens, DdsOverDs = 1.0, center = [0,0]):
    """This is a function to estimate the scale that needs to be used when
plotting the image plane of a gravitational lens. If you're creating a simulated
lens for example, you may not yet know on what scale it should be plotted to
show something useful. This function tries to return a relevant scale that
can be used when specifying the lower-left and upper-right corners for a plot.

Arguments:
 - `lens`: the :py:class:`gravitational lens<grale.lenses.GravitationalLens>`
   that you're interested in.

 - `DdsOverDs`: the fraction Dds/Ds that's relevant for the scenario.

 - `center`: in case the lens is significantly off-center, specifying an estimate
   of the center can help in determining a relevant scale.
"""

    position = ANGLE_DEGREE
    end = ANGLE_ARCSEC/1000.0
    centerX = center[0]
    centerY = center[1]
    
    positions = np.zeros((4,2), dtype = np.double)
    scale = None
    first = True
    factor = 1.2
    
    while position > end:
        positions[0][0] = centerX+position
        positions[0][1] = centerY

        positions[1][0] = centerX-position
        positions[1][1] = centerY

        positions[2][0] = centerX
        positions[2][1] = centerY+position

        positions[3][0] = centerX
        positions[3][1] = centerY-position

        invMag = lens.getInverseMagnification(1.0, DdsOverDs, positions)

        if len(invMag[invMag <= 0]) > 0:
            if first:
                raise PlotException("Can't find scale: starting point already has a negative inverse magnification")
            
            scale = position * factor
            break

        position /= factor
        first = False
    
    if scale is None:
        raise PlotException("Can't find scale: probed positions never appear to have a negative inverse magnification")
    
    plotScaleFactor = 1.5
    scale *= plotScaleFactor
    
    return scale

def _toLensInfo(lensOrLensInfo):
    from . import lenses

    if issubclass(type(lensOrLensInfo), lenses.GravitationalLens):
        return quickLensInfo(lensOrLensInfo)

    if type(lensOrLensInfo) == dict:
        try:
            lensOrLensInfo = DensInfo(**lensOrLensInfo)
        except:
            try:
                lensOrLensInfo = LensInfo(**lensOrLensInfo)
            except Exception as e:
                raise LensInfoException("Can't interpret dictionary as DensInfo or LensInfo")

    return lensOrLensInfo

def quickLensInfo(lens):
    """For the specified `lens`, this creates a `lensInfo` dictionary
that can be used on several plotting functions. The corners of the
region are estimated using :func:`estimatePlotScale`, and `Ds` and
`Dds` entries are set to 1.0.
"""
    s = estimatePlotScale(lens)
    return LensInfo(lens=lens, size = 2.0*s, Dds=1.0, Ds= 1.0)

def calculateDeflectionAndDerivativesForFITS(lens, numXY, angularSize, lensCenterRARec, renderer = "default", 
                                             feedbackObject = "default"):
    """This function calculates a :py:class:`lens plane<grale.images.LensPlane>` for a
specified gravitational lens, and does this in such a way that the
results (deflection angles and their derivatives) can be stored in a 
`FITS <https://en.wikipedia.org/wiki/FITS>`_ file.

Arguments:
 - `lens`: the gravitational lens for which the properties should be calculated.

 - `numXY`: an array of length two containing the number of pixels in X and Y dimentsions
   for the resulting FITS file and grids.

 - `angularSize`: an array of length two containing the angular size in X and Y directions
   respectively.

 - `lensCenterRARec`: for the World Coordinate System in the FITS file, the lens center
   is assumed to be located at these coordinates. This is an array of length two, of which
   the first entry is the right ascension and the second is the declination.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

When successful, the function returns a dictionary containing the following entries:

 - ``emptyfits``: an empty fits file, in which one of the calculated properties could be
   stored.
 - ``alpha_x``: a grid containing the x-component of the deflection angle at the corresponding
   location.
 - ``alpha_xy``: a grid containing the y-component of the deflection angle at the corresponding
   location.
 - ``alpha_xx``: a grid containing the x-derivative of the x-component of the deflection angle.
 - ``alpha_yy``: a grid containing the y-derivative of the y-component of the deflection angle.
 - ``alpha_xy``: a grid containing the y-derivative of the x-component of the deflection angle.
 - ``lensplane``: the lens plane that was calculated. Note that if similar grids would be obtained
   from this lens plane instance directly, they would **not** be the same! Due to the fact that
   the right ascension axis goes from right to left, the values would be mirrored in this
   direction.

Example::

    import grale.lenses as lenses
    import grale.plotutil as plotutil
    import grale.cosmology as cosmology
    from grale.constants import *

    cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)

    # Calculate angular diameter distances for redshift
    zd, zs = 0.5, 1.0
    Dd = cosm.getAngularDiameterDistance(zd)
    Ds = cosm.getAngularDiameterDistance(zs)
    Dds = cosm.getAngularDiameterDistance(zd, zs)

    # Create a lens
    lens = lenses.PlummerLens(Dd, { "width": 2*ANGLE_ARCSEC, "mass": 1e13*MASS_SUN })

    # Create the properties as well as an empty FITS file in which they can be stored
    w = 20*ANGLE_ARCSEC
    ra = 6.6485541367*ANGLE_DEGREE
    dec = 17.1619006374*ANGLE_DEGREE
    r = plotutil.calculateDeflectionAndDerivativesForFITS(lens, [1024, 1024], [ w, w], [ ra, dec], 
                                                          renderer = "openmp")

    # Calculate the convergence (kappa), and write it to a FITS file
    fitsFile = r["emptyfits"]
    axx, ayy = r["alpha_xx"], r["alpha_yy"]

    fitsFile[0].data = 0.5*(Dds/Ds)*(axx+ayy)
    fitsFile.writeto(open("plummer_dens.fits", "wb"))

"""
    numX, numY = numXY
    angularWidth, angularHeight = angularSize
    centerRA, centerDec = lensCenterRARec

    renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "LENSPLANE")

    f = createEmptyFits(numX, numY,
                        1.0*numX/(numX-1)*angularWidth/ANGLE_DEGREE, 
                        1.0*numY/(numY-1)*angularHeight/ANGLE_DEGREE, centerRA/ANGLE_DEGREE,
                        centerDec/ANGLE_DEGREE, numX/2.0, numY/2.0)

    from astropy import wcs

    w = wcs.WCS(f[0].header)
    corner0 = images.centerOnPosition(np.array(w.all_pix2world(0.5,0.5, 0))*ANGLE_DEGREE, [centerRA, centerDec])
    corner1 = images.centerOnPosition(np.array(w.all_pix2world(numX-0.5,numY-0.5, 0))*ANGLE_DEGREE, [centerRA, centerDec])

    diff = corner0-corner1
    #print(diff/ANGLE_ARCSEC)

    bottomLeft = [min(corner0[0], corner1[0]), min(corner0[1], corner1[1])]
    topRight = [max(corner0[0], corner1[0]), max(corner0[1], corner1[1])]

    lp = images.LensPlane(lens, bottomLeft, topRight, numX, numY, renderer = renderer, feedbackObject = feedbackObject)

    alphas = lp.getAlphas()
    ax = np.fliplr(alphas["alpha_x"])
    ay = np.fliplr(alphas["alpha_y"])

    derivs = lp.getAlphaVectorDerivatives()
    axx = np.fliplr(derivs["alpha_xx"])
    ayy = np.fliplr(derivs["alpha_yy"])
    axy = np.fliplr(derivs["alpha_xy"])

    return { 
        "lensplane": lp,
        "emptyfits": f,
        "alpha_x": ax,
        "alpha_y": ay,
        "alpha_xx": axx,
        "alpha_yy": ayy,
        "alpha_xy": axy,
    }

def _getLensFunctionAndDistance(lensInfo):
    lens = lensInfo.getLens()
    Dd = lens.getLensDistance()

    gf = gridfunction.GridFunction(lensInfo.getDensityPoints(), bottomLeft = lensInfo.getBottomLeft(),
                                   topRight = lensInfo.getTopRight())
    return gf.evaluate, Dd

def plotAverageDensityProfile(lensOrLensInfo, thetaMax, center = [0.0, 0.0], thetaSteps = 512, phiSteps = 512, 
                              angularUnit = "default", densityUnit = 1.0, axes = None, renderer = "default",
                              feedbackObject = "default", **kwargs):
    """Creates a plot of the circularly averaged density profile.

Arguments:
 - `lensOrLensInfo`: information about the gravitational lens for which the
   plot should be made. To obtain the lens density values on a grid, for numerical
   integration, the :func:`plotDensity` function is called first (without actually
   plotting the results). There, you can find more information about this
   parameter.

 - `thetaMax`: the angular radius up to which the plot should be made.

 - `center`: the center of the circular averaging.

 - `thetaSteps`: to do the calculations numerically the maximum radius
   `thetaMax` will be divided into this many parts.

 - `phiSteps`: to average the mass density at a specific radius, the
   circular region will be subdivided into this many parts.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `densityUnit`: by default, the density will be in kg/m^2, but another unit can be specified
   here.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. The value `False` has a special meaning: in that case,
   the calculations will be performed as usual, but an actual plot will not be created.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

"""
    angularUnit = _getAngularUnit(angularUnit)

    # Obtain the map if not yet available
    lensInfo = plotDensity(lensOrLensInfo, axes=False, renderer=renderer, feedbackObject=feedbackObject)

    gf = gridfunction.GridFunction(lensInfo.getDensityPoints(), bottomLeft = lensInfo.getBottomLeft(),
                                   topRight = lensInfo.getTopRight())
    F = gf.evaluate
    rMax, rSteps = thetaMax, thetaSteps

    averageValues = np.zeros([rSteps], dtype=np.double)
    r = np.linspace(0, rMax, rSteps)
    pos = np.empty([rSteps, 2], dtype=np.double)

    phiSpace = np.linspace(0, 2*np.pi, phiSteps+1)[:-1] # We don't want both phi = 0 and phi = 2 pi
    for phi in phiSpace:
        pos[:,0] = center[0] + r*np.cos(phi)
        pos[:,1] = center[1] + r*np.sin(phi)
        values = F(pos)
        averageValues += values

    averageValues /= phiSteps

    if axes is not False:
        if axes is None:
            import matplotlib.pyplot as plt
            axes = plt.gca()

        axes.plot(r/angularUnit, averageValues/densityUnit, **kwargs)

    return r, averageValues

def plotIntegratedMassProfile(lensOrLensInfo, thetaMax, center = [0.0, 0.0], thetaSteps = 512, phiSteps = 512,
                              angularUnit = "default", massUnit = 1.0, axes = None, renderer = "default", feedbackObject = "default",
                              **kwargs):
    """Creates a plot of the circularly integrated mass profile. The arguments are 
the same as for :func:`plotAverageDensityProfile`."""

    angularUnit = _getAngularUnit(angularUnit)

    # Obtain the map if not yet available
    lensInfo = plotDensity(lensOrLensInfo, axes=False, renderer=renderer, feedbackObject=feedbackObject)

    F, Dd = _getLensFunctionAndDistance(lensInfo)
    rMax, rSteps = thetaMax, thetaSteps

    integratedValues = np.zeros([rSteps], dtype=np.double)
    r = np.linspace(0, rMax, rSteps+1)
    rFirst = r[:-1]
    rNext = r[1:]
    pos = np.empty([rSteps+1, 2], dtype=np.double)

    for phi in np.linspace(0, 2*np.pi, phiSteps):
        pos[:,0] = center[0] + r*np.cos(phi)
        pos[:,1] = center[1] + r*np.sin(phi)
        values = F(pos)
        vFirst = values[:-1]
        vNext = values[1:]
        integratedValues += rFirst*(2.0*vFirst+vNext) + rNext*(vFirst + 2.0*vNext)

    rLimits = np.linspace(0, rMax, rSteps+1)[1:]
    integratedValues *= 1.0/6.0 * rMax/rSteps * 2.0*np.pi/phiSteps

    integratedValues *= Dd**2

    iv = np.zeros([rSteps])
    iv[0] = integratedValues[0]
    for i in range(1, len(iv)):
        iv[i] = iv[i-1] + integratedValues[i]

    if axes is not False:
        if axes is None:
            import matplotlib.pyplot as plt
            axes = plt.gca()

        axes.plot(rLimits/angularUnit, iv/massUnit, **kwargs)

    return rLimits, iv

def plotSubdivisionGrid(cells, angularUnit = "default", axes = None, **kwargs):
    """Creates a plot of the specified subdivision grid, obtained by
a function from the :mod:`grid<grale.grid>` module for example.

Arguments:
 - `cells`: the grid cells.
 
 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. The value `False` has a special meaning: in that case,
   the calculations will be performed as usual, but an actual plot will not be created.

 - `kwargs`: these parameters will be passed on to the `imshow <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.imshow.html>`_
   function in matplotlib.
"""
    angularUnit = _getAngularUnit(angularUnit)

    from . import grid

    xpoints, ypoints = [], []
    for c in cells["cells"]:
        rc = grid._realCellCenterAndSize(cells["size"], cells["center"], c)
        cx, cy = rc[0]
        w2 = rc[1]/2.0

        cx /= angularUnit
        cy /= angularUnit
        w2 /= angularUnit

        xpoints += [cx-w2, cx-w2, cx+w2, cx+w2, cx-w2, None]
        ypoints += [cy-w2, cy+w2, cy+w2, cy-w2, cy-w2, None]

    if axes is not False:
        if axes is None:
            import matplotlib.pyplot as plt
            axes = plt.gca()

        axes.plot(xpoints, ypoints, '-', **kwargs)

    return [ xpoints, ypoints ]


def plotImagesData(imgDat, angularUnit = "default", plotHull = False, skipTriangles = False, axes = None, fillTriangles = False):
    """Creates a plot of either a single images data set, or a list
of them. In the first case, each image will be drawn with a different
color, in the second, each complete images data set (typically
for a specific source) uses a different color.

Arguments:
 - `imgDat`: either a list, typically to plot several sources at once,
   or a single entry. The entry itself can be an instance
   if :py:class:`ImagesData <grale.images.ImagesData>` or can
   be a dictionary of which the ``imgdata`` key contains such
   an ``ImagesData`` object.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `plotHull`: a flag indicating if the convex hull of the points for
   each image should be plotted as well. Can be useful for visualizing
   which points belong to the same image, in case no triangulation is
   stored in the ``ImagesData`` instance.

 - `skipTriangles`: if this flag is set, triangles stored in the
   ``ImagesData`` instance will not be drawn.

 - `axes`: the default will cause a new plot to be created, but you can specify an existing
   matplotlib axes object as well. 

 - `fillTriangles`: if set, the triangulation, if present, will be filled.
"""
    angularUnit = _getAngularUnit(angularUnit)

    def processImage(im, idx, x, y, segs):
        for p in im.getImagePoints(idx):
            x.append(p["position"][0]/angularUnit)
            x.append(None)
            y.append(p["position"][1]/angularUnit)
            y.append(None)

        if not skipTriangles:
            if im.hasTriangulation():
                for corners in im.getTriangles(idx):
                    corners = list(map(lambda q: tuple(im.getImagePointPosition(idx, q)), corners))

                    seg1 = tuple(sorted([corners[0], corners[1]], key = lambda xy: xy[0]**2 + xy[1]**2))
                    seg2 = tuple(sorted([corners[1], corners[2]], key = lambda xy: xy[0]**2 + xy[1]**2))
                    seg3 = tuple(sorted([corners[2], corners[0]], key = lambda xy: xy[0]**2 + xy[1]**2))

                    segs.add(seg1)
                    segs.add(seg2)
                    segs.add(seg3)

    def getImageTriangles(im, idx, triangles, fillColor):
        if skipTriangles:
            return

        if triangles is None:
            return

        if im.hasTriangulation():
            for corners in im.getTriangles(idx):
                corners = np.array(list(map(lambda q: tuple(im.getImagePointPosition(idx, q)), corners)))/angularUnit

                triangles.append([ corners[0][0], corners[1][0], corners[2][0] ])
                triangles.append([ corners[0][1], corners[1][1], corners[2][1] ])
                triangles.append(fillColor)
            
    def segsToCoords(segs, xsegs, ysegs):

        for s in segs:
            xsegs.append(s[0][0]/angularUnit)
            xsegs.append(s[1][0]/angularUnit)
            xsegs.append(None)
            ysegs.append(s[0][1]/angularUnit)
            ysegs.append(s[1][1]/angularUnit)            
            ysegs.append(None)
        
        return xsegs, ysegs
    
    def addHull(im, idx, x, y):
        
        if im.getNumberOfImagePoints(idx) < 3:
            return
        
        from scipy.spatial import ConvexHull
        xy = np.array([ [ p["position"][0], p["position"][1] ] for p in im.getImagePoints(idx) ])/angularUnit
        hull = ConvexHull(xy)

        for pidx in hull.simplices:
            pt = im.getImagePointPosition(idx, pidx[0])/angularUnit
            x.append(pt[0])
            y.append(pt[1])
            pt = im.getImagePointPosition(idx, pidx[1])/angularUnit
            x.append(pt[0])
            y.append(pt[1])
            x.append(None)
            y.append(None)
    
    import matplotlib.colors as colors
    if axes is None:
        import matplotlib.pyplot as plt
        axes = plt.gca()

    if type(imgDat) == list: # Each source in different color
        for im in imgDat:
            im = im["imgdata"] if type(im) == dict and "imgdata" in im else im
            x, y = [], []
            segs = set()

            num = im.getNumberOfImages()
            for idx in range(num):
                processImage(im, idx, x, y, segs)
                if plotHull:
                    addHull(im, idx, x, y)
            
            segsToCoords(segs, x, y)
            p = axes.plot(x, y, ".-")

            col = p[0].get_color()
            col = colors.to_rgb(col) + (0.3,)
            col = colors.to_hex(col, True)
            triangles = [ ] if fillTriangles else None

            if triangles is not None:
                for idx in range(num):
                    getImageTriangles(im, idx, triangles, col)
                axes.fill(*triangles)
            
    else: # Each image of this one source in different color
        imgDat = imgDat["imgdata"] if type(imgDat) == dict and "imgdata" in imgDat else imgDat
        num = imgDat.getNumberOfImages()
        for idx in range(num):
            x, y = [], []
            segs = set()
            processImage(imgDat, idx, x, y, segs)
                
            segsToCoords(segs, x, y)
            if plotHull:
                addHull(imgDat, idx, x, y)

            p = axes.plot(x, y, ".-")

            col = p[0].get_color()
            col = colors.to_rgb(col) + (0.3,)
            col = colors.to_hex(col, True)
            triangles = [ ] if fillTriangles else None
            if triangles is not None:
                getImageTriangles(imgDat, idx, triangles, col)
                axes.fill(*triangles)

def createEmptyFits(xPix, yPix, angularWidthDeg, angularHeightDeg, raCenterDeg, decCenterDeg, xpCenter, ypCenter):
    """Helper function to create an empty FITS file. 

Arguments:

 - `xPix`: number of pixels in the FITS file in the X-direction
 - `yPix`: number of pixels in the FITS file in the Y-direction
 - `angularWidthDeg`: the angular width of the file, in degrees
 - `angularHeightDeg`: the angular height of the file, in degrees
 - `raCenterDeg`: the right ascension of the `xpCenter` pixel, in degrees
 - `decCenterDeg`: the declination of the center `ypCenter` pixel, in degrees
 - `xpCenter`: the pixel position that corresponds to `raCenterDeg`
 - `ypCenter`: the pixel position that corresponds to `decCenterDeg`
"""
    from astropy.io import fits

    zeros = np.zeros((yPix,xPix),dtype="double")
    hdu = fits.PrimaryHDU(zeros)
    hdu.header["equinox"] = 2000.
    hdu.header["radecsys"] = "FK5     "
    hdu.header["ctype1"]= "RA---TAN"
    hdu.header["crval1"] = raCenterDeg
    hdu.header["crpix1"] = xpCenter+0.5
    hdu.header["cdelt1"] = -1.0*angularWidthDeg/xPix
    hdu.header["cunit1"] = 'deg     '
    hdu.header["ctype2"] = 'DEC--TAN'
    hdu.header["crval2"] = decCenterDeg
    hdu.header["crpix2"] = ypCenter+0.5
    hdu.header["cdelt2"] = 1.0*angularHeightDeg/yPix
    hdu.header["cunit2"] = 'deg     '
    hdu.header["date"] = time.strftime("%Y-%m-%dT%H:%M:%S")

    hdulist = fits.HDUList([hdu])
    return hdulist

def plotImagePlaneFITS(lensInfo, sources, raCenterDeg, decCenterDeg, renderer = "default",
                       feedbackObject = "default", plotSources = False):
    """Creates a FITS image (actually a HDUList that can then be saved) for a
lensing situation.

Arguments:
 - `lensInfo`: an instance of :class:`LensInfo`

 - `sources`: a list of :ref:`source shapes <sourceshapes>` that should be used to calculate
   the images from.

 - `raCenterDeg`: the right ascension of the center of the rendered image plane, in degrees.

 - `decCenterDeg`: the declination of the center of the rendered image plane, in degrees

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

 - `plotSources`: if `True`, the sources themselves will be plotted instead of the images.

"""
    xPix, yPix = lensInfo.getNumXPixels(), lensInfo.getNumYPixels()
    tr, bl = lensInfo.getTopRight(), lensInfo.getBottomLeft()
    angularWidthDeg = (tr[0]-bl[0])/ANGLE_DEGREE
    angularHeightDeg = (tr[1]-bl[1])/ANGLE_DEGREE
    xpCenter = xPix/2.0
    ypCenter = yPix/2.0
    f = createEmptyFits(xPix, yPix, angularWidthDeg, angularHeightDeg, raCenterDeg, decCenterDeg, xpCenter, ypCenter)
    
    imgPlane = lensInfo.getImagePlane(renderer, feedbackObject)
    data = imgPlane.renderImages(sources) if not plotSources else imgPlane.renderSources(sources)
    f[0].data = np.fliplr(data)
    return f
                    
def plotDensityFITS(lensInfo, raCenterDeg, decCenterDeg, renderer = "default",
                    feedbackObject = "default"):
    """Creates a FITS image (actually a HDUList that can then be saved) for the
mass density in `lensInfo`.

Arguments:
 - `lensInfo`: an instance of :class:`LensInfo` or of :class:`DensInfo`.

 - `raCenterDeg`: the right ascension of the center of the rendered image plane, in degrees.

 - `decCenterDeg`: the declination of the center of the rendered image plane, in degrees

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. 

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.
"""
    xPix, yPix = lensInfo.getNumXPixels(), lensInfo.getNumYPixels()
    tr, bl = lensInfo.getTopRight(), lensInfo.getBottomLeft()
    angularWidthDeg = (tr[0]-bl[0])/ANGLE_DEGREE
    angularHeightDeg = (tr[1]-bl[1])/ANGLE_DEGREE
    xpCenter = xPix/2.0
    ypCenter = yPix/2.0
    f = createEmptyFits(xPix, yPix, angularWidthDeg, angularHeightDeg, raCenterDeg, decCenterDeg, xpCenter, ypCenter)
    
    data = lensInfo.getDensityPixels(renderer, feedbackObject)
    f[0].data = np.fliplr(data)
    return f

def _getAngularUnit(u):
    if u == "default":
        return getDefaultAngularUnit()
    return u

_defaultAngularUnit = [ 1.0 ]

def getDefaultAngularUnit():
    """In various plot functions you can specify which `angularUnit` should be
    used when visualizing the results. If it is set to ``"default"``, then the
    value returned by this function will be used. It can be set with
    :func:`setDefaultAngularUnit`."""
    return _defaultAngularUnit[0]

def setDefaultAngularUnit(x):
    """In various plot functions you can specify which `angularUnit` should be
    used when visualizing the results. If it is set to ``"default"``, then the
    value set by this function will be used."""
    _defaultAngularUnit[0] = x

