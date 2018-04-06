from __future__ import print_function

def hoursMinutesSecondsToDegrees(s):
    """Converts a string or array of three parts, specifying hours, minutes and seconds,
    into a single floating point number that corresponds to a number of degrees.

    As an example, passing ``"01:23:45"``, ``"01 23 45"``, ``["01","23","45"]`` and
    ``[1,23,45]`` will all produce the same output of 20.9375.
    """

    if type(s) != list:
        if ":" in s:
            s = s.split(":")
        else: # Perhaps a space is used?
            s = s.split(" ")
            
    if len(s) != 3:
        raise Exception("Expected three parts, but detected {}".format(len(s)))
        
    h, m, s = float(s[0]), float(s[1]), float(s[2])
    
    if h < 0 or h >= 24:
        raise Exception("'Hours' must lie between 0 and 24, but is {}".format(h))
    if m < 0 or m >= 60:
        raise Exception("'Minutes must lie between 0 and 60, but is {}".format(m))
    if s < 0 or s >= 60:
        raise Exception("'Seconds must lie between 0 and 60, but is {}".format(s))
    
    return (((s / 60.0) + m)/60.0 + h)/24.0 * 360.0

def degreesMinutesSecondsToDegrees(s):
    """Converts a string or array of three parts, specifying degrees, minutes and seconds,
    into a single floating point number that corresponds to a number of degrees.

    As an example, passing ``"-1:23:45"``, ``"-1 23 45"``, ``["-1","23","45"]`` and
    ``[-1,23,45]`` will all produce the same output of -1.39583.
    """
    if type(s) != list:
        if ":" in s:
            s = s.split(":")
        else: # Perhaps a space is used
            s = s.split(" ")
            
    if len(s) != 3:
        raise Exception("Expected three parts, but detected {}".format(len(s)))
        
    d, m, s = float(s[0]), float(s[1]), float(s[2])
    
    sign = 1.0
    if d < 0:
        sign = -1.0
        d = -d

    if m < 0 or m >= 60:
        raise Exception("'Minutes must lie between 0 and 60, but is {}".format(m))
    if s < 0 or s >= 60:
        raise Exception("'Seconds must lie between 0 and 60, but is {}".format(s))
    
    return (((s / 60.0) + m)/60.0 + d) * sign

# For internal use
def _processImagesDataDict(imgs):
    if not imgs:
        raise Exception("Empty images set")

    from .images import ImagesData
    
    imgDat = ImagesData(len(imgs))
    groupInfo = { }
    
    imgIdx = 0
    for imgId in imgs:
        imgPoints = imgs[imgId]
        if len(imgPoints) == 0:
            raise Exception("Empty image")
            
        for ptInfo in imgPoints:
            ptIdx = imgDat.addPoint(imgIdx, [ ptInfo["x"], ptInfo["y"] ])
            if "group" in ptInfo:
                grpId = ptInfo["group"]
                if not grpId in groupInfo:
                    groupInfo[grpId] = []
                groupInfo[grpId].append((imgIdx, ptIdx))
                
            if "timedelay" in ptInfo:
                imgDat.addTimeDelayInfo(imgIdx, ptIdx, ptInfo["timedelay"])
        
        imgIdx += 1
        
    # Process the stored groups
    for grpKey in groupInfo:
        grpIdx = imgDat.addGroup()
        for imgIdx, ptIdx in groupInfo[grpKey]:
            imgDat.addGroupPoint(grpIdx, imgIdx, ptIdx)
    
    return imgDat


def _defLA(l):
    from . import constants
    p = l.split()
    info = { 
        "x": float(p[0])*constants.ANGLE_ARCSEC,
        "y": float(p[1])*constants.ANGLE_ARCSEC
    }
    if len(p) > 2:
        info["z"] = float(p[2])
    return info

_defaultLineAnalyzer = [ _defLA ]

def getDefaultLineAnalyzer():
    """Returns the current default line analyzer, used when the `lineAnalyzer`
    parameter in :func:`readInputImagesFile` is set to ``"default"``."""
    return _defaultLineAnalyzer[0]

def setDefaultLineAnalyzer(x):
    """Sets the current default line analyzer to ``x``, which will be used when the `lineAnalyzer`
    parameter in :func:`readInputImagesFile` is set to ``"default"``."""
    _defaultLineAnalyzer[0] = x

def _getLineAnalyzer(la):
    if la == "default":
        return getDefaultLineAnalyzer()
    return la

def readInputImagesFile(inputData, isPointImagesFile, lineAnalyzer = "default", centerOn = [0, 0]):
    """This function can process a text file (or previously read text data) into
    one or more :class:`ImagesData` instances.

    If the function you pass in ``lineAnalyzer`` returns a source identifier
    and image identifier, these will be used to group images per source. Otherwise,
    the behaviour is different depending on the value of ``isPointImagesFile``. If
    ``True``, then each single line is interpreted as a single point image, and a
    blank line groups the point images per source. If ``isPointImagesFile`` is ``False``,
    then a single blank line separates the points that belong to different images,
    and a double blank line separates the images from different sources.

    Arguments:
     - `inputData`: this can be an open file object, a filename or just
       the text data that's already read from a file. Lines that start
       with a ``#`` sign are considered to be comments, and are ignored
       completely.

     - `isPointImagesFile`: as explained above, this flag indicates how input
       is treated when no source and image identifiers are specified in the
       `lineAnalyzer` function.

     - `lineAnalyzer`: here you can pass a function that interprets a single
       (non-empty) line in the input data. As default, the first column will
       be interpreted as the x-coordinate, the second as the y-coordinate,
       both expressed in arcseconds. The third column, if present will be used
       as the redshift. A different default can be set using
       :func:`setDefaultLineAnalyzer`.
       
       If you provide the function yourself, it should return a dictionary with
       these entries:
       
        - ``x``: the x-coordinate of the point, converted to radians (the
          :ref:`pre-defined constants <constants>` can be useful here).
          The helper functions :func:`degreesMinutesSecondsToDegrees` and
          :func:`hoursMinutesSecondsToDegrees` can also be of assistance.
        - ``y``: the y-coordinate of the point, converted to radians (the
          :ref:`pre-defined constants <constants>` can be useful here)
          The helper function :func:`degreesMinutesSecondsToDegrees` can also be of assistance.
        - (optionally) ``srcnr``: the source identifier of this point.
        - (optionally) ``imgnr``: the image identifier of this point.
        - (optionally) ``z``: the redshift for the source. Note that if this
          is specified for more than one image in a source, they should all
          have exactly the same redshift.
        - (optionally): ``group``: if certain points in different images
          belong together, this can be specified by having the same group
          identifier.
        - (optionally): ``timedelay``: if timedelay information is known for
          a point, it can be passed along this way.

     - `centerOn`: if present, all coordinates will be recentered on this
       point, using the :func:`ImagesData.centerOnPosition` function.

    The return value of this function is a list of dictionaries with entries:

     - ``imgdata``: the :class:`ImagesData` object for a specific source
     - (optionally) ``z``: if it was present in the input, the redshift for
       the source is stored.
     - (optionally) ``srcnr``: if a specific source identifier was specified
       in the file, and passed along in the `lineAnalyzer` function, it is
       stored here.

    Examples: suppose that each line contains entries like the one
    below (from `A free-form lensing grid solution for A1689 with new multiple images <http://adsabs.harvard.edu/abs/2015MNRAS.446..683D>`_) ::

        # i  ID   B05  REF  RAJ2000(h:m:s)  DECJ2000(d:m:s)  z     Delta_beta
        1    1.1  1.1  B05  13:11:26.257    -1:19:58.753     3.04  1.03
        2    1.2  1.2  B05  13:11:26.088    -1:20:02.261     3.04  0.73
        3    1.3  1.3  B05  13:11:29.584    -1:21:09.475     3.04  2.50
        ...

    Then, you could use the following function as the `lineAnalyzer` parameter:

    .. code-block:: python

        def la(l):
            p = l.split()
            if len(p) != 8:
                raise Exception("Expecting 8 input fields in line '{}'".format(l))
            
            info = { }
            info["srcnr"], info["imgnr"] = map(int, p[1].split("."))
            info["z"] = p[6]
            info["x"] = hoursMinutesSecondsToDegrees(p[4])*ANGLE_DEGREE
            info["y"] = degreesMinutesSecondsToDegrees(p[5])*ANGLE_DEGREE
            
            return info

    If we'd like to center the coordinates on that first point, we could calculate

    .. code-block:: python

        centerRa = hoursMinutesSecondsToDegrees("13:11:26.257")*ANGLE_DEGREE
        centerDec = degreesMinutesSecondsToDegrees("-1:19:58.753")*ANGLE_DEGREE

    and pass ``[centerRa, centerDec]`` as the `centerOn` parameter.

    As another example, in `Dark matter dynamics in Abell 3827: new data consistent with standard Cold Dark Matter <http://adsabs.harvard.edu/abs/2017arXiv170804245M>`_
    you'll find coordinates like the ones below::

        # Name  RA         Dec
        Ao1     330.47479  -59.94358
        Ao2     330.46649  −59.94665
        ...


    The 'A' is a label for the source, and as all points refer to the same image, is present
    for each point. The second character, 'o' in this case, refers to a feature in an image.
    The third part, '1' and '2' above, indicates the image the point is a part of. We can
    treat this input in two ways: the image number goes from 1 to 7, so we can process the
    input in the following way to create seven images, where each point is part of a specific
    group:

    .. code-block:: python

        def la(l):
            p = l.split()
            if len(p) != 3:
                raise Exception("Expecting 3 input fields in '{}'".format(l))
            
            info = { }
            info["srcnr"] = 'A'
            info["imgnr"] = p[0][2]
            info["group"] = p[0][1]
            info["x"] = float(p[1])*ANGLE_DEGREE
            info["y"] = float(p[2])*ANGLE_DEGREE
            
            return info

    Alternatively, we can treat each feature, each set of corresponding points (marked with 'o'
    for example), as a different point source. This is wat the following `lineAnalyzer`
    function would do:

    .. code-block:: python

        def la(l):
            p = l.split()
            if len(p) != 3:
                raise Exception("Expecting 3 input fields in '{}'".format(l))
            
            info = { }
            info["srcnr"] = p[0][1]
            info["imgnr"] = p[0][2]
            info["x"] = float(p[1])*ANGLE_DEGREE
            info["y"] = float(p[2])*ANGLE_DEGREE
            
            return info
    """

    lineAnalyzer = _getLineAnalyzer(lineAnalyzer)

    try:
        # Check if it's a file that's already open
        lines = inputData.readlines()
    except:
        try:
            with open(inputData, "rt") as f:
                lines = f.readlines()
        except:
            # Check if it's a string that can be split in lines
            lines = inputData.splitlines()

    internalIdPrefix = "justsomelongandunlikelyprefix_"
    internalSourceId = 0
    internalImageId = 0
    sourceCount = 0 # To be able to keep the same ordering of the sources as in the file
    
    def idStr(x):
        return "{}_{:05d}".format(internalIdPrefix, x)
    
    emptyLineCount = 0
    sourceData = { }
    
    for l in lines:
        l = l.strip()
        if l.startswith("#"): # a commented line, ignore completely
            continue
        
        if not l:
            emptyLineCount += 1
        else:
            emptyLineCount = 0
          
        if emptyLineCount == 0: # New info
            if isPointImagesFile:
                internalImageId += 1
            
            info = lineAnalyzer(l)
            imagesId = idStr(internalImageId) if not "imgnr" in info else info["imgnr"]
            sourceId = idStr(internalSourceId) if not "srcnr" in info else info["srcnr"]
            
            if not sourceId in sourceData:
                imagesData = { "points": { } }
                sourceData[sourceId] = imagesData
                sourceData[sourceId]["sourcecount"] = sourceCount
                sourceCount += 1
            else:
                imagesData = sourceData[sourceId]
                
            if not imagesId in imagesData["points"]:
                imagePoints = [ ]
                imagesData["points"][imagesId] = imagePoints
            else:
                imagePoints = imagesData["points"][imagesId]
            
            if "z" in info:
                if "z" in imagesData:
                    if info["z"] != imagesData["z"]:
                        raise Exception("Error in line '{}': redshift '{}' is not the same as previous recorded redshift '{}'".format(info["z"],imagesData["z"]))
                
                imagesData["z"] = info["z"]
            
            ptInfo = {
                "x": float(info["x"]),
                "y": float(info["y"])
            }
            if "group" in info:
                ptInfo["group"] = info["group"]
            if "timedelay" in info:
                ptInfo["timedelay"] = info["timedelay"]
            
            imagePoints.append(ptInfo)
    
        elif emptyLineCount == 1 and not isPointImagesFile:
            internalImageId += 1
        elif emptyLineCount == 2 or (emptyLineCount == 1 and isPointImagesFile):
            internalSourceId += 1
            internalImageId = 0
    
    srcListTmp = sorted([ (sourceData[s]["sourcecount"], s) for s in sourceData])
        
    srcList = [ ]
    for tmpCount, srcId in srcListTmp:
        info = { }
        if "z" in sourceData[srcId]:
            info["z"] = sourceData[srcId]["z"]
        if not str(srcId).startswith(internalIdPrefix):
            info["srcnr"] = srcId
        
        info["imgdata"] = _processImagesDataDict(sourceData[srcId]["points"])
        info["imgdata"].centerOnPosition(centerOn[0], centerOn[1])
        srcList.append(info)
    
    # TODO: if point images file, check that every images data set is indeed a point image set
    #       if not a point image set, _optionally_ check that they're not point images
    
    #return sourceData
    return srcList

def enlargePolygon(points, offset, simplifyScale = 0.02):
    """TODO"""
    from shapely.geometry.polygon import LinearRing
    from .images import ImagesDataException

    if len(points) < 4:
        raise ImagesDataException("At least four points must be present")
    if points[0][0] != points[-1][0] or points[0][1] != points[-1][1]:
        raise ImagesDataException("Not closed: first point should equal last point")

    def getScale(x):
        minx, miny, maxx, maxy = x.bounds
        scale = ((maxx-minx)**2 + (maxy-miny)**2)**0.5
        return scale

    r = LinearRing(points[:-1])
    if offset < 0: # Negative indicates fractional enlarge
        offset = (-offset)*getScale(r)
        
    direction = "right" if r.is_ccw else "left"
    o = r.parallel_offset(offset, direction, join_style=2)
    if hasattr(o, "geoms"):
        raise Exception("Couldn't find a single shape after adding border")
    
    o = LinearRing(o.coords)
    simpl = o.simplify(simplifyScale*getScale(o), preserve_topology=False) if simplifyScale > 0 else o
    pts = [ pt for pt in simpl.coords ]
    return pts + pts[:1]

def createGridTriangles(bottomLeft, topRight, numX, numY, holes = None, enlargeHoleOffset = None,
                        simplifyScale = 0.02, triangleExe = "triangle"):
    """TODO"""
    import tempfile
    import os
    import numpy as np
    import subprocess
    import copy
    import matplotlib.path as mplPath
    from .images import ImagesData, ImagesDataException

    holes = [] if not holes else copy.deepcopy(holes)
    
    filesToDelete = [ ]
    f = tempfile.NamedTemporaryFile("w+t", suffix=".poly", delete=False)
    filesToDelete.append(f.name)
    eleName = f.name[:-4] + "1.ele"
    
    try:
        X,Y = np.meshgrid(np.linspace(bottomLeft[0], topRight[0], numX), 
                          np.linspace(bottomLeft[1], topRight[1], numY))

        XY = np.empty(X.shape + (2,))
        XY[:,:,0], XY[:,:,1]= X, Y
        XY = XY.reshape((-1,2))

        for hIdx in range(len(holes)):
            if enlargeHoleOffset is None:
                raise ImagesDataException("When holes are present, you always need to specify the amount with which they should be enlarged")

            h = holes[hIdx]
            if enlargeHoleOffset != 0:
                holes[hIdx] = enlargePolygon(h, enlargeHoleOffset, simplifyScale)

            if len(h) < 4:
                raise ImagesDataException("Hole must contain at least four points")

            #print(h[0])
            if h[0][0] != h[-1][0] or h[0][1] != h[-1][1]:
                raise ImagesDataException("Hole is not closed!")

        holesPoints = sum([ len(h)-1 for h in holes])
        totalPoints = holesPoints + len(XY)
        pointIds = { }

        ptIdx = 1      
        polyData = "{} 2 0 0\n".format(totalPoints)

        for hIdx in range(len(holes)):
            h = holes[hIdx]
            hNew = [ ]

            for i in range(len(h)-1):
                pointIds[ptIdx] = { "xy": h[i] }
                polyData += "    {} {:.15g} {:.15g}\n".format(ptIdx, h[i][0], h[i][1])
                hNew.append({ "idx": ptIdx, "xy": h[i] })
                ptIdx += 1

            hNew.append(hNew[0]) # Close it
            holes[hIdx] = hNew

        for i in range(len(XY)):
            pointIds[ptIdx] = { "xy": XY[i] }
            polyData += "    {} {:.15g} {:.15g}\n".format(ptIdx, XY[i,0], XY[i,1])
            ptIdx += 1

        # Lines of the holes
        polyData += "{} 0\n".format(holesPoints)
        lnIdx = 1
        for h in holes:
            for i in range(len(h)-1):
                polyData += "    {} {} {}\n".format(lnIdx, h[i]["idx"], h[i+1]["idx"])
                lnIdx += 1

        def getInternalHolePoint(h):
            bbPath = mplPath.Path(np.array([d["xy"] for d in h]))
            for i in range(len(h)-1):
                x0,y0 = h[i]["xy"] 
                x2, y2 = h[(i+2)%(len(h)-1)]["xy"]
                x, y = (x0+x2)/2.0, (y0+y2)/2.0
                if bbPath.contains_point([x,y]):
                    return x, y

            raise ImagesDataException("No internal point found")

        # Internal points for the holes
        polyData += "{}\n".format(len(holes))
        hIdx = 1
        for h in holes:
            x, y = getInternalHolePoint(h)
            polyData += "{} {:.15g} {:.15g}\n".format(hIdx, x, y)
            hIdx += 1

        #print(polyData)
        f.write(polyData)
        f.close()

        if "GRALE_DEBUG_TRIANGLE" in os.environ:
            DEVNULL = None
            print(polyData)
        else:
            try:
                DEVNULL = subprocess.DEVNULL
            except AttributeError:
                DEVNULL = open(os.devnull, "wb")
        
        subprocess.check_call([triangleExe, "-pcNP", f.name ], stdout=DEVNULL, stderr=subprocess.STDOUT)
        filesToDelete.append(eleName)
            
        triangData = open(eleName, "rt").read().splitlines()
        numTriangles, dummy1, dummy2 = list(map(int, triangData[0].split()))
        #print(triangData)

        img = ImagesData(1)
        for i in range(numTriangles):
            tIdx, pt1, pt2, pt3 = list(map(int, triangData[i+1].split()))
            imgIdx = [ ]

            for pt in [ pt1, pt2, pt3 ]:
                e = pointIds[pt]
                if not "imgidx" in e:
                    e["imgidx"] = img.addPoint(0, e["xy"])

                imgIdx.append(e["imgidx"])

            img.addTriangle(0, imgIdx[0], imgIdx[1], imgIdx[2])

        return img
    finally:
        try:
            f.close()
        except Exception as e:
            print("Warning: couldn't close f: {}".format(e))
            
        for n in filesToDelete:
            try:
                os.unlink(n)
            except Exception as e:
                print("Warning: couldn't remove '{}': {}".format(n, e)) 

