from PyQt5 import QtWidgets, QtCore, QtGui
import copy
from grale.constants import * # TODO
from base import Layer, PointGraphicsItemBase, LayerGraphicsItemBase
from debug import log


class MultiplePointGraphicsItem(LayerGraphicsItemBase):
    def __init__(self, layer, isPointSelect = False, parent = None):
        ptType = LayerGraphicsItemBase.PointSelect if isPointSelect else LayerGraphicsItemBase.Normal
        super(MultiplePointGraphicsItem, self).__init__(layer, ptType, False, parent)

class PointsLayer(Layer):
    def __init__(self, name = "Points layer", isPointSelect = False):
        self.isPointSelect = isPointSelect
        super(PointsLayer, self).__init__(name)

    def importFromImagesData(self, imgDat, imgIndex, triangles = True, timedelays = True, groups = True):
        allPointsAdded, allTrianglesAdded = [], []
        if imgIndex == "all":
            for imgIndex in range(imgDat.getNumberOfImages()):
                ptsAdded, triangAdded = self._importFromImagesDataSingle(imgDat, imgIndex, triangles, timedelays, groups)
                allPointsAdded += ptsAdded
                allTrianglesAdded += triangAdded
        else:
            allPointsAdded, allTrianglesAdded = self._importFromImagesDataSingle(imgDat, imgIndex, triangles, timedelays, groups)
        return allPointsAdded, allTrianglesAdded

    def _importFromImagesDataSingle(self, imgDat, imgIndex, triangles, timedelays, groups):
        pts = [ { "xy": imgDat.getImagePointPosition(imgIndex, pt)/ANGLE_ARCSEC,
                  "uuid": self.generateKey() } 
                for pt in range(imgDat.getNumberOfImagePoints(imgIndex)) ]
        
        pointsAdded = set()
        trianglesAdded = set()

        if timedelays:
            for i in range(imgDat.getNumberOfTimeDelays()):
                imgNr, pointNr, delay = imgDat.getTimeDelay(i)
                if imgNr == imgIndex:
                    pts[pointNr]["timedelay"] = delay

        if groups:
            for g in range(imgDat.getNumberOfGroups()):
                for p in range(imgDat.getNumberOfGroupPoints(g)):
                    imgNr, pointNr = imgDat.getGroupPointIndices(g, p)
                    if imgNr == imgIndex:
                        pts[pointNr]["label"] = str(g)

        if triangles and imgDat.hasTriangulation():
            triangList = { }
            for p0, p1, p2 in imgDat.getTriangles(imgIndex):
                triangKey = self.generateKey()
                triangList[triangKey] = set([pts[p]["uuid"] for p in [p0, p1, p2]])

        for p in range(len(pts)):
            key, xy = pts[p]["uuid"], pts[p]["xy"]
            td = None if "timedelay" not in pts[p] else pts[p]["timedelay"]
            label = None if "label" not in pts[p] else pts[p]["label"]

            self.setPoint(key, [float(xy[0]), float(xy[1])], label, td)
            pointsAdded.add(key)

        if triangles and imgDat.hasTriangulation():
            # Keep track of the triangles a point's involved in. This helps
            # when deleting a point
            for triangKey in triangList:
                self.addTriangle(*triangList[triangKey], triangKey)
                trianglesAdded.add(triangKey)
        
        return pointsAdded, trianglesAdded

    def createGraphicsItem(self):
        return MultiplePointGraphicsItem(self, self.isPointSelect)

    def toSettings(self):
        pts = self.getPoints()
        ptIdx = 0
        ptKeyIdx = { }
        for k in pts:
            ptKeyIdx[k] = ptIdx
            ptIdx += 1

        triangles = self.getTriangles();

        settings = {
                "name": self.getName(),
                "points": [ pts[k] for k in pts ],
                "triangles": [ list(map(lambda x: ptKeyIdx[x], triangles[t])) for t in triangles ]
        }

        return settings

    @staticmethod
    def fromSettings(settings):
        l = PointsLayer(settings["name"])
        keys = [ ]
        for p in settings["points"]:
            keys.append(l.addPoint(p["xy"], p["label"], p["timedelay"]))

        for pts in settings["triangles"]:
            pts = list(map(lambda x: keys[x], pts))
            l.addTriangle(*pts)
            
        return l
