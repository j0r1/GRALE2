"""This module defines functions to plot the mass density of a gravitational
lens, or the lens effect itself. The module also contains classes to help
you make animations.
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

class PlotException(Exception):
    """An exception that will be thrown in case something goes wrong when plotting"""
    pass

def _prepareDensity(lensInfo, renderer, feedbackObject):

    renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "MASSDENS")

    if "numx" not in lensInfo:
        lensInfo["numx"] = 511
    if "numy" not in lensInfo:
        lensInfo["numy"] = 511

    bottomLeft = lensInfo["bottomleft"]
    topRight = lensInfo["topright"]
    numX = lensInfo["numx"] + 1
    numY = lensInfo["numy"] + 1

    if "densitypixels" in lensInfo and "densitypoints" in lensInfo:

        massMap = lensInfo["densitypoints"]
        pixels = lensInfo["densitypixels"]

        if numX == massMap.shape[1] and numY == massMap.shape[0]:

            if numX == pixels.shape[1]+1 and numY == pixels.shape[0]+1:
                feedbackObject.onStatus("Reusing density pixels and points")
                return
            else:
                lensInfo["densitypixels"] = None
        else:
            lensInfo["densitypoints"] = None

    if "densitypoints" not in lensInfo:

        lens = lensInfo["lens"]
        massMap = lens.getSurfaceMassDensityMap(bottomLeft, topRight, numX, numY, renderer = renderer, feedbackObject = feedbackObject, reduceToPixels = False)
        pixels = 0.25*massMap[0:numY-1,0:numX-1] + 0.25*massMap[0:numY-1,1:numX] + 0.25*massMap[1:numY,0:numX-1] + 0.25*massMap[1:numY,1:numX]

        lensInfo["densitypixels"] = pixels
        lensInfo["densitypoints"] = massMap

    else: # densitypoints in lensInfo, but not densitypixels

        massMap = lensInfo["densitypoints"]
        pixels = 0.25*massMap[0:numY-1,0:numX-1] + 0.25*massMap[0:numY-1,1:numX] + 0.25*massMap[1:numY,0:numX-1] + 0.25*massMap[1:numY,1:numX]
        lensInfo["densitypixels"] = pixels

    massPixelSum = sum(sum(lensInfo["densitypixels"]))
    totalWidth = abs(topRight[0] - bottomLeft[0])
    totalHeight = abs(topRight[1] - bottomLeft[1])

    Dd = None
    if "Dd" in lensInfo:
        Dd = lensInfo["Dd"]
    elif "lens" in lensInfo:
        lens = lensInfo["lens"]
        Dd = lens.getLensDistance()
        lensInfo["Dd"] = Dd

    if Dd is not None:
        lensInfo["totalmass"] = Dd**2 * (totalWidth/lensInfo["numx"])*(totalHeight/lensInfo["numy"]) * massPixelSum

    return lensInfo


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

     - `lensOrLensInfo`: this can either be a gravitational lens instance or dictionary 
       that contains at least the following entries

        - `lens`: the gravitational lens that should be plotted
        - `bottomleft`: the bottom-left corner of the plot region
        - `topright`: the top-right corner of the plot region

       optionally, the number of pixels that should be plotted can be specified as well:

        - `numx`: number of pixels in x-direction, defaults to 511
        - `numy`: number of pixels in y-direction, defaults to 511

       If possible, the calculated mass density will be used to calculate the total mass
       within the region, which will be stored in a `totalmass` entry of the dictionary.
       The dictionary will also be used to cache some calculated properties in.

       In case it's only a gravitational lens, an estimate of the corners will be used.

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
       to speed up the calculation. If left to ``None``, the default, single core renderer 
       will be used.

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
    _prepareDensity(lensInfo, renderer, feedbackObject)

    # TODO: can we replace this using a gridfunction?
    plotArray = privutilcython.resample2DArray(lensInfo["densitypoints"], numY, numX)
    bottomLeft = lensInfo["bottomleft"]
    topRight = lensInfo["topright"]

    pixels = lensInfo["densitypixels"]
    bottomLeft = lensInfo["bottomleft"]
    topRight = lensInfo["topright"]

    X, Y = np.meshgrid(np.linspace(bottomLeft[0], topRight[0], numX),
                       np.linspace(bottomLeft[1], topRight[1], numY))

    X /= angularUnit
    Y /= angularUnit
    Z = plotArray/densityUnit

    
    plot3DInteractive(X, Y, Z, flipX=flipX, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, initialCamera=initialCamera, 
                      visJSoptions=visJSoptions, maxPoints=maxPoints, height=height, visJS=visJS, visCSS=visCSS,
                      canvas2svgJS=canvas2svgJS, fileSaverJS=fileSaverJS)
    return lensInfo

def plotDensity(lensOrLensInfo, renderer = "default", feedbackObject = "default", angularUnit = "default", densityUnit = 1.0,
                axes = None, axImgCallback = None, **kwargs):
    """Creates a 2D plot of the situation specified in `lensOrLensInfo`.

In case you just want the calculations to be performed (for example because you need
the calculated density pixels) but you don't want an actual plot, you can set the
`axes` parameter to `False`.

Arguments:

 - `lensOrLensInfo`: this can either be a gravitational lens instance or dictionary 
   that contains at least the following entries

    - `lens`: the gravitational lens that should be plotted
    - `bottomleft`: the bottom-left corner of the plot region
    - `topright`: the top-right corner of the plot region

   optionally, the number of pixels that should be plotted can be specified as well:

    - `numx`: number of pixels in x-direction, defaults to 511
    - `numy`: number of pixels in y-direction, defaults to 511

   If possible, the calculated mass density will be used to calculate the total mass
   within the region, which will be stored in a `totalmass` entry of the dictionary.
   The dictionary will also be used to cache some calculated properties in.

   In case it's only a gravitational lens, an estimate of the corners will be used.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. If left to `None`, the default, single core renderer 
   will be used.

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
    _prepareDensity(lensInfo, renderer, feedbackObject)

    if axes is not False:

        pixels = lensInfo["densitypixels"]
        bottomLeft = lensInfo["bottomleft"]
        topRight = lensInfo["topright"]

        if axes is None:
            import matplotlib.pyplot as plt
            axes = plt.gca()

        # Note: need to swap Y labeling here because of the way the pixels are ordered in this function
        axImg = axes.imshow(pixels/densityUnit, extent = np.array([ bottomLeft[0], topRight[0], topRight[1], bottomLeft[1]])/angularUnit, **kwargs)
        axes.invert_yaxis()
        if axImgCallback: axImgCallback(axImg)

    return lensInfo

def _prepareImagePlane(lensInfo, renderer, feedbackObject, evenError, cosmology):

    bottomLeft = lensInfo["bottomleft"]
    topRight = lensInfo["topright"]
    
    if "numx" not in lensInfo:
        lensInfo["numx"] = 511
    if "numy" not in lensInfo:
        lensInfo["numy"] = 511

    if evenError and (lensInfo["numx"] % 2 == 0 or lensInfo["numy"] % 2 == 0):
        raise PlotException("Using an odd number of pixels may work better for simple lenses, set evenError = False to disable this check")

    renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "LENSPLANE")

    imgPlane = None
    if "imageplane" in lensInfo:
        imgPlane = lensInfo["imageplane"]
        if type(imgPlane) == images.ImagePlane:
            Ds = lensInfo["Ds"]
            Dds = lensInfo["Dds"]

            if imgPlane.getDs() != Ds or imgPlane.getDds() != Dds:
                imgPlane = None
        else:
            # TODO: check against lenses?
            pass

        if imgPlane: # Make sure the render region and resolution hasn't changed
            renderInfo = imgPlane.getRenderInfo()
            if renderInfo["bottomleft"][0] != bottomLeft[0] or renderInfo["bottomleft"][1] != bottomLeft[1] or \
               renderInfo["topright"][0] != topRight[0] or renderInfo["topright"][1] != topRight[1] or \
               renderInfo["xpixels"] != lensInfo["numx"] or renderInfo["ypixels"] != lensInfo["numy"]:
                imgPlane = None

    if imgPlane is None:
        lensPlane = None
        if "lensplane" in lensInfo:
            lensPlane = lensInfo["lensplane"]
            renderInfo = lensPlane.getRenderInfo()

            if renderInfo["bottomleft"][0] != bottomLeft[0] or renderInfo["bottomleft"][1] != bottomLeft[1] or \
               renderInfo["topright"][0] != topRight[0] or renderInfo["topright"][1] != topRight[1] or \
               renderInfo["xpoints"] != lensInfo["numx"]+1 or renderInfo["ypoints"] != lensInfo["numy"]+1:
                lensPlane = None

        if lensPlane is None:
            if type(lensInfo["lens"]) != list:
                if cosmology is not None:
                    raise PlotException("For a single lens, the 'cosmology' parameter is not used and should be set to 'None'")

                lensPlane = images.LensPlane(lensInfo["lens"], bottomLeft, topRight, lensInfo["numx"]+1, lensInfo["numy"]+1, renderer = renderer, feedbackObject = feedbackObject)
            else: # multiple lens
                #if renderer is not None:
                #    raise PlotException("No special renderers are currently supported for a multi-lensplane scenario, the renderer should be set to 'None'")
                if cosmology is None:
                    raise PlotException("For a multiple lensplane scenario, the cosmology parameter must be set")

                lensPlane = multiplane.MultiLensPlane(lensInfo["lens"], cosmology, bottomLeft, topRight, lensInfo["numx"]+1, lensInfo["numy"]+1, renderer = renderer, feedbackObject = feedbackObject)

            lensInfo["lensplane"] = lensPlane
        else:
            feedbackObject.onStatus("Reusing lens plane")

        if type(lensInfo["lens"]) != list:
            Ds = lensInfo["Ds"]
            Dds = lensInfo["Dds"]
            imgPlane = images.ImagePlane(lensPlane, Ds, Dds)
        else:
            zs = lensInfo["zs"]
            imgPlane = multiplane.MultiImagePlane(lensPlane, zs)

        lensInfo["imageplane"] = imgPlane
    else:
        feedbackObject.onStatus("Reusing image plane")

    return lensInfo
    
def plotImagePlane(lensOrLensInfo, sources = [], renderer = "default", feedbackObject = "default", 
                   angularUnit = "default", subSamples = 9, sourceRgb = (0, 1, 0), imageRgb = (1, 1, 1),
                   plotCaustics = True, plotCriticalLines = True, plotSources = True, plotImages = True,
                   evenError = True, axes = None, cosmology = None, axImgCallback = None, **kwargs):
    """Create a matplotlib-based plot of image plane and/or source plane for certain
lens parameters in `lensOrLensInfo`. You can also use this function to create the necessary
calculated mappings but not the plot, by setting `axes` to `False`. This can be useful
to create a :py:class:`LensPlane<grale.images.LensPlane>` and :py:class:`ImagePlane<grale.images.ImagePlane>` instance
using a specific renderer to speed up the calculation.

Arguments:

 - `lensOrLensInfo`: this can either be a gravitational lens instance or dictionary 
   that contains at least the following entries

    - `lens`: the gravitational lens that used for the plot
    - `bottomleft`: the bottom-left corner of the plot region
    - `topright`: the top-right corner of the plot region
    - `Ds`: the angular diameter distance to the source
    - `Dds`: the angular diameter distance between lens and source

   optionally, the number of pixels that should be plotted can be specified as well:

    - `numx`: number of pixels in x-direction, defaults to 511
    - `numy`: number of pixels in y-direction, defaults to 511

   Upon completion, this dictionary will contain `lensplane` and `imageplane`, containing
   a :py:class:`LensPlane<grale.images.LensPlane>` and :py:class:`ImagePlane<grale.images.ImagePlane>` instance 
   respectively.

   In case it's only a gravitational lens, an estimate of the corners will be used, and
   both `Ds` and `Dds` will be set to 1.0 (only the value of Dds/Ds matters, and will be
   one in this case).

 - `sources`: a list of :ref:`source shapes <sourceshapes>` that should be used to calculate
   the images from.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. If left to ``None``, the default, single core renderer 
   will be used.

 - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.

 - `angularUnit`: the angular unit that should be used in the plot. The 
   :ref:`pre-defined constants <constants>` can be useful here.

 - `subSamples`: each pixel is sub-sampled sqrt(subSamples) times in x- and y- direction, to be able to 
   roughly approximate the integration that’s needed over the surface area of a pixel. This parameter
   is passed to :py:meth:`ImagePlane.renderImages<grale.images.ImagePlane.renderImages>`
   and :py:meth:`ImagePlane.renderSources<grale.images.ImagePlane.renderSources>`.

 - `sourceRgb`: the RGB color (each component is a number between 0 and 1) that should be
   given to a pixel that lies within the source.

 - `imageRgb`: the RGB color (each component is a number between 0 and 1) that should be
   given to a pixel that lies within an image.

 - `plotCaustics`: boolean value that indicates if the caustics should be drawn on the plot.

 - `plotCriticalLines`: boolean value that indicates if the critical lines should be drawn on
   the plot.

 - `plotSources`: boolean value that indicates if the sources should be drawn on the plot.

 - `plotImages`: boolean value that indicates if the images should be drawn on the plot.

 - `evenError`: by default, the function will raise an exception if the number of pixels
   in the x or y-direction is even, because for simple lenses odd pixel numbers work better.
   In case you don't care for this, set the flag to `False` and this will stop the calculation.

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
    _prepareImagePlane(lensInfo, renderer, feedbackObject, evenError, cosmology)

    bottomLeft = lensInfo["bottomleft"]
    topRight = lensInfo["topright"]
    imgPlane = lensInfo["imageplane"]
    numX = lensInfo["numx"]
    numY = lensInfo["numy"]

    if axes is None:
        import matplotlib.pyplot as plt
        axes = plt.gca()

    # TODO: sourceScale and imageScale
    
    rgbSplane = np.zeros((numY, numX, 3))
    rgbIplane = np.zeros((numY, numX, 3))

    if sources:

        if plotSources:
            splane = imgPlane.renderSources(sources, subSamples = subSamples)

            rgbSplane[:,:,0] = splane*sourceRgb[0]
            rgbSplane[:,:,1] = splane*sourceRgb[1]
            rgbSplane[:,:,2] = splane*sourceRgb[2]

        if plotImages:
            iplane = imgPlane.renderImages(sources, subSamples = subSamples)
            rgbIplane[:,:,0] = iplane*imageRgb[0]
            rgbIplane[:,:,1] = iplane*imageRgb[1]
            rgbIplane[:,:,2] = iplane*imageRgb[2]

    if axes is not False:
        # Note: need to swap Y labeling here because of the way the pixels are ordered in this function
        axImg = axes.imshow(rgbIplane + rgbSplane, extent = np.array([ bottomLeft[0], topRight[0], topRight[1], bottomLeft[1]])/angularUnit, **kwargs)
        if axImgCallback: axImgCallback(axImg)
        
    criticalLines = imgPlane.getCriticalLines()

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

    if plotCaustics:
        caustics = imgPlane.getCaustics()
        plotLines(caustics, angularUnit, color="blue")

    if plotCriticalLines:
        criticalLines = imgPlane.getCriticalLines()
        plotLines(criticalLines, angularUnit, color="red")

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
                          flipX = False, gnuplotExe = "gnuplot", convertExe = "convert", cosmology = None):

    """Creates a gnuplot-based plot of the image plane/source plane that corresponds to the
situation specified in `lensOrLensInfo`, saving various files with names starting with 
`fileNameBase`. If you only want to generate these files, but not visualize them in a plot 
immediately, you can set `axes` to `False`.

Arguments:

 - `lensOrLensInfo`: this can either be a gravitational lens instance or dictionary 
   that contains at least the following entries

    - `lens`: the gravitational lens that used for the plot
    - `bottomleft`: the bottom-left corner of the plot region
    - `topright`: the top-right corner of the plot region
    - `Ds`: the angular diameter distance to the source
    - `Dds`: the angular diameter distance between lens and source

   optionally, the number of pixels that should be plotted can be specified as well:

    - `numx`: number of pixels in x-direction, defaults to 511
    - `numy`: number of pixels in y-direction, defaults to 511

   Upon completion, this dictionary will contain `lensplane` and `imageplane`, containing
   a :py:class:`LensPlane<grale.images.LensPlane>` and :py:class:`ImagePlane<grale.images.ImagePlane>` instance 
   respectively.

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
   to speed up the calculation. If left to ``None``, the default, single core renderer 
   will be used.

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

    _prepareImagePlane(lensInfo, renderer, feedbackObject, evenError, cosmology)

    bottomLeft = lensInfo["bottomleft"]
    topRight = lensInfo["topright"]
    imgPlane = lensInfo["imageplane"]
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

 - `lensOrLensInfo`: this can either be a gravitational lens instance or dictionary 
   that contains at least the following entries

    - `lens`: the gravitational lens that should be plotted
    - `bottomleft`: the bottom-left corner of the plot region
    - `topright`: the top-right corner of the plot region

   optionally, the number of pixels that should be plotted can be specified as well:

    - `numx`: number of pixels in x-direction, defaults to 511
    - `numy`: number of pixels in y-direction, defaults to 511

   If possible, the calculated mass density will be used to calculate the total mass
   within the region, which will be stored in a `totalmass` entry of the dictionary.
   The dictionary will also be used to cache some calculated properties in.

   In case it's only a gravitational lens, an estimate of the corners will be used.

 - `fileNameBase`: the generated gnuplot file will have this name, ending with ``.gnuplot``.
   Similarly, the generated eps file will start with this and end with ``.eps``. If set to
   `None`, no output files will be generated, but the result may still be shown.

 - `numX`: the number of points in the x-direction used to create the surface or 2D
   plot.

 - `numY`: the number of points in the y-direction used to create the surface or 2D
   plot.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. If left to ``None``, the default, single core renderer 
   will be used.

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
    _prepareDensity(lensInfo, renderer, feedbackObject)
    
    # TODO: can we replace this using a gridfunction?
    plotArray = privutilcython.resample2DArray(lensInfo["densitypoints"], numY, numX)
    bottomLeft = lensInfo["bottomleft"]
    topRight = lensInfo["topright"]

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
    return lensOrLensInfo

def quickLensInfo(lens):
    """For the specified `lens`, this creates a `lensInfo` dictionary
that can be used on several plotting functions. The corners of the
region are estimated using :func:`estimatePlotScale`, and `Ds` and
`Dds` entries are set to 1.0.
"""
    s = estimatePlotScale(lens)
    return { "lens": lens, "bottomleft": [ -s, -s], "topright": [s, s], 
             "Dd": lens.getLensDistance(), "Dds": 1.0, "Ds": 1.0 }

def calculateDeflectionAndDerivativesForFITS(lens, numXY, angularSize, lensCenterRARec, renderer = "default", 
                                             feedbackObject = "default"):
    """This function calculates a :py:class:`lens plane<grale.images.LensPlane>` for a
specified gravitational lens, and does this in such a way that the
results (deflection angles and their derivatives) can be stored in a 
`FITS <https://en.wikipedia.org/wiki/FITS>`_ file.

Arguments:

 - ``lens``: the gravitational lens for which the properties should be calculated.

 - ``numXY``: an array of length two containing the number of pixels in X and Y dimentsions
   for the resulting FITS file and grids.

 - ``angularSize``: an array of length two containing the angular size in X and Y directions
   respectively.

 - ``lensCenterRARec``: for the World Coordinate System in the FITS file, the lens center
   is assumed to be located at these coordinates. This is an array of length two, of which
   the first entry is the right ascension and the second is the declination.

 - `renderer`: this parameter can be used to specify a specific :ref:`renderer <renderers>`
   to speed up the calculation. If left to ``None``, the default, single core renderer 
   will be used.

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

    def createEmptyFits(xPix, yPix, angularWidthDeg, angularHeightDeg, raCenterDeg, decCenterDeg, xpCenter, ypCenter):

        from astropy.io import fits

        zeros = np.zeros((yPix,xPix),dtype="double")
        hdu = fits.PrimaryHDU(zeros)
        hdu.header["equinox"] = 2000.
        hdu.header["radecsys"] = "FK5     "
        hdu.header["ctype1"]= "RA---TAN"
        hdu.header["crval1"] = raCenterDeg
        hdu.header["crpix1"] = xpCenter+1
        hdu.header["cdelt1"] = -1.0*angularWidthDeg/xPix
        hdu.header["cunit1"] = 'deg     '
        hdu.header["ctype2"] = 'DEC--TAN'
        hdu.header["crval2"] = decCenterDeg
        hdu.header["crpix2"] = ypCenter+1
        hdu.header["cdelt2"] = 1.0*angularHeightDeg/yPix
        hdu.header["cunit2"] = 'deg     '
        hdu.header["date"] = time.strftime("%Y-%m-%dT%H:%M:%S")

        hdulist = fits.HDUList([hdu])
        return hdulist

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

def _getLensFunctionAndDistance(lensOrLensInfo):
    #if type(lensOrLensInfo) == dict:

    lens = lensOrLensInfo["lens"]
    Dd = lens.getLensDistance()

    gf = gridfunction.GridFunction(lensOrLensInfo["densitypoints"], bottomLeft = lensOrLensInfo["bottomleft"],
                                   topRight = lensOrLensInfo["topright"])
    F = gf.evaluate

    #else:
    #    lens = lensOrLensInfo
    #    F = lens.getSurfaceMassDensity
    #    Dd = lens.getLensDistance()

    return F, Dd

def plotAverageDensityProfile(lensOrLensInfo, thetaMax, center = [0.0, 0.0], thetaSteps = 512, phiSteps = 512, 
                              angularUnit = "default", densityUnit = 1.0, axes = None, renderer = "default",
                              feedbackObject = "default", **kwargs):
    """Creates a plot of the circularly averaged density profile

TODO
"""
    angularUnit = _getAngularUnit(angularUnit)

    # Obtain the map if not yet available
    lensInfo = plotDensity(lensOrLensInfo, axes=False, renderer=renderer, feedbackObject=feedbackObject)

    F, Dd = _getLensFunctionAndDistance(lensInfo)
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
    """Creates a plot of the circularly integrated mass profile.

TODO
"""
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
    """Creates a plot of the specified subdivision grid

TODO
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


def plotImagesData(imgDat, angularUnit = "default", plotHull = False, skipTriangles = False, axes = None):
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
            axes.plot(x, y, ".-")
            
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

            axes.plot(x, y, ".-")
                    
def _getAngularUnit(u):
    if u == "default":
        return getDefaultAngularUnit()
    return u

_defaultAngularUnit = [ 1.0 ]

def getDefaultAngularUnit():
    return _defaultAngularUnit[0]

def setDefaultAngularUnit(x):
    _defaultAngularUnit[0] = x
