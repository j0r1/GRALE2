import copy
import pprint
from pointslayer import PointsLayer

import grale.images as images # TODO
from grale.constants import ANGLE_ARCSEC
import cppqt

def strToBool(s):
    if type(s) == bool:
        return s
    if s.lower() == "false":
        return False
    if s.lower() == "true":
        return True
    raise Exception("Unrecognized boolean string '{}'".format(s))

def valueFromSettings(settings, name, castFunction, defaultValue):
    x = settings.value(name)
    x = castFunction(x) if x is not None else defaultValue
    return x

# Now implemented in C++, to speed things up for a null space for example
###def _splitPointsAndTriangles(allPoints, allTriangles, pointsLeft = None):
###
###    imageInfo = [ ]
###
###    trianglesInPoints = { }
###    for t in allTriangles:
###        for p in allTriangles[t]:
###            if not p in trianglesInPoints:
###                trianglesInPoints[p] = set()
###            trianglesInPoints[p].add(t)
###
###    while True:
###        curPoints = { }
###        curTriangles = { }
###
###        t = None
###        for k in allTriangles:
###            t = k
###            break
###
###        if not t: # Can't find a triangle to start from
###            break
###
###        # Start with this triangle, and these points
###        curTriangles[t] = allTriangles[t]
###        del allTriangles[t]
###
###        #print("Found triangle", t)
###        #pprint.pprint(curTriangles[t])
###
###        for p in curTriangles[t]:
###            curPoints[p] = allPoints[p]
###            del allPoints[p]
###
###        # Look for other triangles in these points:
###        while True:
###            
###            t = None
###            for p in curPoints:
###                for k in trianglesInPoints[p]:
###                    # Check if we've already processed this triangle 
###                    if k in allTriangles:
###                        t = k
###                        break
###
###                if t:
###                    break
###
###            if not t: # No remaining triangle found, stop
###                break
###
###            # Otherwise add the triangle and it's points
###            curTriangles[t] = allTriangles[t]
###            del allTriangles[t]
###            #print("Found next triangle", t)
###            #pprint.pprint(curTriangles[t])
###
###            for p in curTriangles[t]:
###                if p in allPoints: # It's possible that the point was already present in another triangle
###                    curPoints[p] = allPoints[p]
###                    del allPoints[p]
###
###        imageInfo.append({
###            "points": curPoints,
###            "triangles": curTriangles
###        })
###
###    # Check if there are still points left
###    for p in allPoints:
###        if pointsLeft is not None:
###            for p2 in allPoints:
###                pointsLeft.append(p2)
###        raise Exception("Can't split points layer in multiple images: still points left after considering all triangulations")
###
###    return imageInfo

def layersToImagesData(layers, multipleImagesPerLayer = True, saveGroups = True, saveTimeDelays = True,
                       pointsLeftInfo = None):
    
    imageInfo = [ ]

    for l in layers:
        if type(l) != PointsLayer:
            continue

        if not multipleImagesPerLayer:
            imageInfo.append({
                "points": l.getPoints(),
                "triangles": l.getTriangles()
            })
        else:
            if pointsLeftInfo is not None:
                d = { "layer": l.getUuid(), "pointsleft": [] }
                pointsLeftInfo.append(d)
            #imageInfo += _splitPointsAndTriangles(l.getPoints(), l.getTriangles(), pointsLeft=pl)

            r = cppqt.cppqt.splitPointsAndTriangles(l)
            pprint.pprint(r)
            if type(r) == str: # an error message
                raise Exception(r)
            if type(r) == dict:
                d["pointsleft"] = r["pointsleft"]
                raise Exception(r["error"])

            imageInfo += r

    if len(imageInfo) == 0:
        raise Exception("No image information found that can be exported to an images data object")

    # First add these points and triangles to an ImagesData instance,
    # keeping track of the image and point id
    imgDat = images.ImagesData(len(imageInfo))

    imagePoints = [ { } for i in range(len(imageInfo)) ]

    for imgIdx in range(imgDat.getNumberOfImages()):
        points = imageInfo[imgIdx]["points"]

        for p in points:
            pt = points[p]
            #pprint.pprint(pt)
            ptIdx = imgDat.addPoint(imgIdx, ( pt["xy"][0]*ANGLE_ARCSEC, pt["xy"][1]*ANGLE_ARCSEC ))
            imagePoints[imgIdx][ptIdx] = pt
            pt["indices"] = (imgIdx, ptIdx)
    
    # Add the triangles
    for imgIdx in range(imgDat.getNumberOfImages()):
        points = imageInfo[imgIdx]["points"]
        triangles = imageInfo[imgIdx]["triangles"]
        
        for t in triangles:
            ptIndices = [ ]
            for pKey in triangles[t]:
                imgIdx2, ptIdx = points[pKey]["indices"]
                assert(imgIdx2 == imgIdx) # Make sure it's all from the same image

                ptIndices.append(ptIdx)

            assert(len(set(ptIndices)) == 3) # Make sure it's three different points

            imgDat.addTriangle(imgIdx, *ptIndices) 

    # Check that point group labels are unique per image, and build group info
    groupInfo = { }
    for points in imagePoints:
        groups = set()
        for ptIdx in points:
            pt = points[ptIdx]
            groupName = pt["label"]
            if not groupName:
                continue

            if groupName in groups:
                raise Exception("Same group name '{}' is present more than once within a single image".format(groupName))

            if not groupName in groupInfo:
                groupInfo[groupName] = [ ]
            groupInfo[groupName].append(pt["indices"])

    # Store the groups
    if saveGroups:
        for groupName in groupInfo:
            gId = imgDat.addGroup()
            points = groupInfo[groupName]
            for imIdx, ptIdx in points:
                imgDat.addGroupPoint(gId, imIdx, ptIdx)

    # Store the time delay info
    if saveTimeDelays:
        for points in imagePoints:
            for ptIdx in points:
                pt = points[ptIdx]
                if pt["timedelay"] is not None:
                    imIdx, ptIdx = pt["indices"]
                    imgDat.addTimeDelayInfo(imIdx, ptIdx, pt["timedelay"])

    return imgDat

def main():
    #imgDat = images.ImagesData.load("/tmp/multimages.imgdata")
    #layer = PointsLayer()
    #layer.importFromImagesData(imgDat, 0)
    import json
    layer = PointsLayer.fromSettings(json.load(open("test.json")))

    newImgDat = layersToImagesData([layer])

    import grale.plotutil as plotutil
    import matplotlib.pyplot as plt
    plotutil.plotImagesData(newImgDat)
    plt.show()

    newImgDat.save("/tmp/splittest.imgdat")

if __name__ == "__main__":
    main()
