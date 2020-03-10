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

def _getLinesFromInputData(inputData):
    try:
        # Check if it's a file that's already open
        return inputData.readlines()
    except:
        pass

    try:
        with open(inputData, "rt") as f:
            return f.readlines()
    except:
        pass

    # Check if it's a string that can be split in lines
    return inputData.splitlines()

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
    lines = _getLinesFromInputData(inputData)

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
                        raise Exception("Error in line '{}': redshift '{}' is not the same as previous recorded redshift '{}'".format(l, info["z"], imagesData["z"]))
                
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

# NOTE: for join_style = 2 the function below can cause some odd spikes in the
#       enlarged polygon, e.g. with:
#        import grale.images as images
#        import numpy as np
#        import matplotlib.pyplot as plt
#        
#        pts = np.array([
#            [-14.468749999999998, 8.2265625],
#            [-14.515625, 8.171875],
#            [-14.59375, 8.140625],
#            [-14.5, 8.1171875],
#            [-14.875, 8.046875],
#            [-14.9921875, 8.046875],
#            [-15.0703125, 8.0859375],
#            [-15.117187499999998, 8.125],
#            [-15.132812499999998, 8.1796875],
#            [-15.093750000000002, 8.3046875],
#            [-14.9921875, 8.359375],
#            [-14.7109375, 8.390625],
#            [-14.476562499999998, 8.359375],
#            [-14.445312500000002, 8.2734375],
#            [-14.468749999999998, 8.2265625]])
#        
#        pts2 = np.array(images.enlargePolygon(pts, offset=0.5))
#        
#        plt.plot(pts[:,0], pts[:,1])
#        plt.plot(pts2[:,0], pts2[:,1])
#        plt.show()

def enlargePolygon(points, offset, simplifyScale = 0.02):
    """For a given set of points that describe a polygon, enlarge this
    by adding a specified offset. The `Shapely <https://shapely.readthedocs.io/en/latest/>`_
    library is used to accomplish this.
    
    Arguments:
     * `points`: this list of points describes the polygon; the last point 
       should be the same as the first point.
     * `offset`: if this value is positive, it describes a distance that's 
       added to the sides of the polygon; if it is negative the absolte 
       value is interpreted as a fraction and the distance that's added is 
       calculated as this fraction times the scale of the polygon. 
     * `simplifyScale`: if greate than zero, the newly obtained polygon is 
       simplified, and this describes a tolerance below which points can be
       removed. It is specified as a fraction of the scale of the polygon.
    """
    for i in range(10, -1, -1):
        try:
            return _enlargePolygon(points, offset, simplifyScale)
        except Exception as e:
            if i == 0:
                raise
            print(f"WARNING({i}): {e}")

    raise Exception("Unexpected: shouldn't get here")

def _enlargePolygon(points, offset, simplifyScale = 0.02):
    from shapely.geometry.polygon import LinearRing, Polygon
    import shapely.ops
    from .images import ImagesDataException

    if points[0][0] != points[-1][0] or points[0][1] != points[-1][1]:
        raise ImagesDataException("Not closed: first point should equal last point")

    def getScale(x):
        minx, miny, maxx, maxy = x.bounds
        scale = ((maxx-minx)**2 + (maxy-miny)**2)**0.5
        return scale

    import copy
    import numpy as np
    points = copy.deepcopy(points)

    for j in range(10,-1,-1): # Try at most 10 times
        points1 = points[:-1]
        r = LinearRing(points1)
        eps = getScale(r)*1e-5

        points2 = points1[len(points1)//2:] + points1[:len(points1)//2]
        points2 = copy.deepcopy(points2)

        for i in range(len(points2)):
            points1[i][0] += np.random.normal(scale=eps)
            points1[i][1] += np.random.normal(scale=eps)
            points2[i][0] += np.random.normal(scale=eps)
            points2[i][1] += np.random.normal(scale=eps)

        r = LinearRing(points1)
        r2 = LinearRing(points2)

        del points1
        del points2
        
        if offset < 0: # Negative indicates fractional enlarge
            offset = (-offset)*(getScale(r)+getScale(r2))*0.5

        try:
            o, o2 = [ x.parallel_offset(offset, "right" if x.is_ccw else "left", join_style=1) for x in [r, r2] ]
            break
        except Exception as e:
            if j == 0:
                raise

            print(f"WARNING({j}): {e}")

    unionObjs = []
    def addPoly(g):
        pts = [ [g.xy[0][i], g.xy[1][i]] for i in range(len(g.xy[0])) ]
        if len(pts) == 2: # add a point in between to create a nearly flat triangle
            cx, cy = (g.xy[0][0] + g.xy[0][1])*0.5, (g.xy[1][0] + g.xy[1][1])*0.5
            vx, vy = g.xy[0][0] - g.xy[0][1], g.xy[1][0] - g.xy[1][1]
            f = 0.001 # TODO: what's a good value here?

            px, py = cx - f*vy, cy + f*vx
            pts.append( [ px, py ])

        # NOTE: to work around a TopologyException here, the buffer(0)
        #       trick is used. See https://stackoverflow.com/a/20873812/2828217
        unionObjs.append(Polygon(pts).buffer(0))

    for x in o, o2:
        if hasattr(x, "geoms"):
            for g in x:
                addPoly(g)
        else:
            addPoly(x)

    o = shapely.ops.unary_union(unionObjs)
    o = LinearRing(o.boundary.coords)

#    if hasattr(o, "geoms"):
#        #print("O:", o)
#        for g in o.geoms:
#            import matplotlib.pyplot as plt
#            plt.plot(g.xy[0], g.xy[1], '.-')
#
#        h = o.convex_hull
#        print(h.exterior)
#        #print("convex hull:", h)
#        plt.plot(h.xy[0], h.xy[1])
#
#        plt.plot(r.xy[0], r.xy[1], 'o-')
#        plt.show()
#
#        raise Exception("Couldn't find a single shape after adding border")
    
    simpl = o.simplify(simplifyScale*getScale(o), preserve_topology=False) if simplifyScale > 0 else o
    pts = [ pt for pt in simpl.coords ]
    return pts + pts[:1]

def _shapelyPlot(o, label=None, show=False):
    import matplotlib.pyplot as plt
    
    if type(o) == list:
        for x in o:
            _shapelyPlot(x)
    else:
        if hasattr(o, "geoms"):
            raise Exception("TODO")
        else:
            plt.plot(o.xy[0], o.xy[1])

    if label:
        plt.gca().set_title(label)

    if show:
        plt.show()

def _checkHoleOverlap(holes):
    from shapely.geometry.polygon import LinearRing, Polygon
    import shapely.ops
    import copy

    holes = [ LinearRing(h[:-1]) for h in holes ]

    #_shapelyPlot(holes, "Holes", show=True)

    # Check if some holes are completely interior to another
    insideSomething = set()
    for i in range(0, len(holes)-1):
        for j in range(i+1, len(holes)):
            if holes[i].within(holes[j]):
                insideSomething.add(i)
            elif holes[j].within(holes[i]):
                insideSomething.add(j)

    # Remove the shapes that are completely inside another
    holes = [ holes[i] for i in range(len(holes)) if not i in insideSomething ]

    # Check intersections
#    intersections = { }
#    for i in range(0, len(holes)-1):
#        for j in range(i+1, len(holes)):
#            if holes[i].intersects(holes[j]):
#                s = set() if not i in intersections else intersections[i]
#                s.add(i)
#                s.add(j)
#                origSet = copy.copy(s)
#                for k in origSet:
#                    if k in intersections:
#                        s |= intersections[k]
#
#                for k in origSet:
#                    intersections[k] = s
    
    
    intersections = { i:set([i]) for i in range(len(holes)) }
    for i in range(0, len(holes)-1):
        for j in range(i+1, len(holes)):
            if holes[i].intersects(holes[j]):
                intersections[i].add(j)
                intersections[j].add(i)

    changed = True
    while changed:
        #import pprint
        #print("Intersections:")
        #pprint.pprint(intersections)
        changed = False

        for i in [ x for x in intersections]:
            origSet = intersections[i].copy()
            for j in origSet:
                intersections[i] |= intersections[j]
                if origSet != intersections[i]:
                    changed = True

    #import pprint
    #print("Final Intersections:")
    #pprint.pprint(intersections)

    #_shapelyPlot(holes, "Intersections", show=True)

    # First add the holes that have no intersections
    newHoles = [ holes[i] for i in range(len(holes)) if not i in intersections ]

    #pprint.pprint(holes)

    # If there are intersections, try to do something sensible
    # For all the things that intersect, we're going to add a union
    for i in intersections:
        if not holes[i]: # We've already processed this
            continue

        # join the intersections in this list
        s = intersections[i]

        #print("Handling i = ", i)
        # Try to use the union function from, and use convex hull if it fails

        totalObj = shapely.ops.unary_union([ Polygon(holes[j].coords) for j in s if holes[j]])
        try:
            totalObj = LinearRing(totalObj.boundary.coords)
        except Exception as e:
            print("Warning: couldn't find union ({}), falling back to convex hull".format(e))
            
            #allPoints = [ pt for pt in holes[j].coords for j in s ]
            totalObj = LinearRing(totalObj.convex_hull.boundary.coords)


        # Mark holes as processed
        for j in s:
            holes[j] = None

        newHoles.append(totalObj)
        #pprint.pprint(totalObj)

    #_shapelyPlot(newHoles, "New holes", True)

    # Convert back to list of points
    holes = [ ]
    for h in newHoles:
        pts = [ pt for pt in h.coords ]
        holes.append(pts + pts[:1])

    return holes

def _getHolesList(holes, gridScale):
    from .images import ImagesData, ImagesDataException
    import numpy as np
    import copy

    newHoles = []
    if holes:
        for h in holes:
            # For an ImagesData instance, try to use the triangulated data
            # to get a boundary for the region, or if not available, use 
            # the convex hull of the points
            if type(h) == ImagesData:
                for i in range(h.getNumberOfImages()):
                    foundBorder = False
                    try:
                        b = h.getBorder(i)
                        foundBorder = True
                    except ImagesDataException as e:
                        #print("Warning: ignoring exception: {}".format(e))
                        pass

                    if foundBorder:
                        newHoles.append(b)
                        continue

                    try:
                        b = h.getConvexHull(i)
                        foundBorder = True
                    except ImagesDataException as e:
                        #print("Warning: ignoring exception: {}".format(e))
                        pass

                    if foundBorder:
                        newHoles.append(b)
                        continue

                    # Couldn't get convex hull or border, should be only one or two points
                    pts = [ p["position"] for p in h.getImagePoints(i)]
                    numPoints = 16 # TODO: configurable?
                    if len(pts) == 1:
                        center = pts[0]
                        A = gridScale/1000 # TODO: configurable?
                        B = A
                        startAngle = 0
                        angles = np.linspace(0, 2*np.pi, numPoints+1)[:-1]
                    elif len(pts) == 2:
                        center = (pts[0] + pts[1])/2
                        A = np.sum((pts[0]-center)**2)**0.5
                        B = A/20 # TODO: configurable?
                        diff = pts[0]-center
                        startAngle = np.arctan2(diff[1], diff[0])
                        angles = np.linspace(0, 2*np.pi, numPoints+1)[:-1]
                    else:
                        raise ImagesDataException(f"Unable to get image border based on either triangulation or convex hull, and number of points ({len(pts)}) is different than expected")

                    coords = np.empty((numPoints, 2))
                    coords[:,0] = (A*np.cos(angles) * np.cos(startAngle) - B*np.sin(angles) * np.sin(startAngle)) + center[0]
                    coords[:,1] = (A*np.cos(angles) * np.sin(startAngle) + B*np.sin(angles) * np.cos(startAngle)) + center[1]
                    b = [xy for xy in coords]
                    newHoles.append(b + b[:1])

            else: # If not an ImagesData instance, assume it's polygon coords
                newHoles.append(copy.deepcopy(h))

    return newHoles

def createGridTriangles(bottomLeft, topRight, numX, numY, holes = None, enlargeHoleOffset = None,
                        simplifyScale = 0.02, triangleExe = "triangle",
                        checkOverlap = True):
    """Creates a grid of triangles, out of which some holes may be cut. This grid can
    then be used as a null space grid in lens inversions. When such holes are cut out,
    it is usually a good idea to make them somewhat larger than the images themselves.
    Even if no enlargement of the holes is required, it must still be specified. 

    The `Triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ program is used
    to create the triangulation.

    The result is returned as an :class:`ImagesData <grale.images.ImagesData>`
    instance.

    Arguments:
     * `bottomLeft`: bottom-left corner of the triangulated region.
     * `topRight`: top-right corner of the triangulated region.
     * `numX`: number of points in the X-direction.
     * `numY`: number of points in the Y-direction.
     * `holes`: a list of holes that should be cut out of the triangulation. 
       If an entry of this list is an :class:`ImagesData <grale.images.ImagesData>
       instance, then for each image a hole will be added, preferrably based on the
       triangulated data that may be stored in the instance (using
       :func:`ImagesData.getBorder <grale.images.ImagesData.getBorder>`), or as a 
       fall-back based on the convex hull of the image (using
       :func:`ImagesData.getConvexHull <grale.images.ImagesData.getConvexHull>`).
       If it is not an ``ImagesData`` instance, it is assumed to be a list of 
       points describing a polygon, where the first point in the
       list must equal the last point.
     * `enlargeHoleOffset`: specifies the value by which the holes need to be
       enlarged. For each hole, this is passed as the `offset` argument of the
       :func:`enlargeHoleOffset` function. If holes are present, it _must_ be
       specified, so if you do not want to enlarge the holes you should set this
       to 0.
     * `simplifyScale`: this is passed to the :func:`enlargeHoleOffset` function
       which is used to enlarge the specified holes.
     * `triangleExe`: the executable for the `Triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ 
       program
    """
    import tempfile
    import os
    import numpy as np
    import subprocess
    import copy
    import matplotlib.path as mplPath
    from .images import ImagesData, ImagesDataException

    gridScale = ((bottomLeft[0]-topRight[0])**2 + (bottomLeft[1]-topRight[1])**2)**0.5
    holes = _getHolesList(holes, gridScale)
    
    filesToDelete = [ ]
    f = tempfile.NamedTemporaryFile("w+t", suffix=".poly", delete=False)
    filesToDelete.append(f.name)
    eleName = f.name[:-4] + "1.ele"
    nodeName = f.name[:-4] + "1.node"
    
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

        if checkOverlap:
            holes = _checkHoleOverlap(holes)

        holesPoints = sum([ len(h)-1 for h in holes])
        totalPoints = holesPoints + len(XY)

        ptIdx = 1      
        polyData = "{} 2 0 0\n".format(totalPoints)

        for hIdx in range(len(holes)):
            h = holes[hIdx]
            hNew = [ ]

            for i in range(len(h)-1):
                polyData += "    {} {:.15g} {:.15g}\n".format(ptIdx, h[i][0], h[i][1])
                hNew.append({ "idx": ptIdx, "xy": h[i] })
                ptIdx += 1

            hNew.append(hNew[0]) # Close it
            holes[hIdx] = hNew

        for i in range(len(XY)):
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
        
        subprocess.check_call([triangleExe, "-pcP", f.name ], stdout=DEVNULL, stderr=subprocess.STDOUT)
        filesToDelete.append(eleName)
        filesToDelete.append(nodeName)
            
        # Read the points that are used in the triangulation file
        pointIds = { }
        nodeData = open(nodeName, "rt").read().splitlines()
        numPts, dummy1, dummy2, dummy3 = list(map(int, nodeData[0].split()))
        if "GRALE_DEBUG_TRIANGLE" in os.environ:
            print("\n".join(nodeData))

        for i in range(numPts):
            ptId, x, y, dummy = list(map(float, nodeData[i+1].split()))
            pointIds[ptId] = { "xy": [ x, y] }

        # Read the triangles
        triangData = open(eleName, "rt").read().splitlines()
        numTriangles, dummy1, dummy2 = list(map(int, triangData[0].split()))
        if "GRALE_DEBUG_TRIANGLE" in os.environ:
            print("\n".join(triangData))

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
            
        if not "GRALE_DEBUG_TRIANGLE_KEEPFILES" in os.environ:
            for n in filesToDelete:
                try:
                    os.unlink(n)
                except Exception as e:
                    print("Warning: couldn't remove '{}': {}".format(n, e)) 

def createSourceFromImagesData(imgDat, idx = -1):
    """For the points in the :class:`ImagesData <grale.images.ImagesData>` instance,
    a :class:`source shape <grale.images.SourceImage>` will be created. This is
    useful after obtaining backprojected images using e.g. 
    :func:`InversionWorkSpace.backProject <grale.inversion.InversionWorkSpace.backProject>`,
    to get an estimate of the source shape, and re-calculate the image positions.
    If only one point is available, a :class:`PointSource <grale.images.PointSource>`
    instance will be created, otherwise a :class:`PolygonSource <grale.images.PolygonSource>`
    is used. If `idx` is negative, all points are used, otherwise a specific image
    is selected.
    """
    if idx < 0:
        points = [ imgDat.getImagePointPosition(i,j) for i in range(imgDat.getNumberOfImages()) for j in range(imgDat.getNumberOfImagePoints(i)) ]
    else:
        points = [ imgDat.getImagePointPosition(idx,j) for j in range(imgDat.getNumberOfImagePoints(idx)) ]

    if not points:
        raise Exception("No points were specified in the images data instance")

    from .images import PolygonSource, PointSource
    import numpy as np

    if len(points) == 1:
        return PointSource(points[0])

    if len(points) == 2:
        # Add some points to make a small diamond shape
        cx = (points[0][0]+points[1][0])*0.5
        cy = (points[0][1]+points[1][1])*0.5
        vx = points[0][0] - cx
        vy = points[0][1] - cy
        f = 0.1 # TODO: what's a good value here?
        wx = -vy*f
        wy = vx*f
        pt1 = [ cx+wx, cy+wy ]
        pt3 = [ cx-wx, cy-wy ]
        points = points[:1] + [ pt1 ] + points[1:] + [ pt3 ]

    points = np.array(points)
    position = points[0,:]
    points = points - position
    return PolygonSource(position, points, True)

def createPointImagesData(thetas):
    """This creates an :class:`ImagesData <grale.images.ImagesData>` instance,
    with one point per image, according to the entries in the `thetas` list
    of positions."""
    from .images import ImagesData
    imgDat = ImagesData(len(thetas))
    for tIdx in range(len(thetas)):
        imgDat.addPoint(tIdx, thetas[tIdx])
    return imgDat
