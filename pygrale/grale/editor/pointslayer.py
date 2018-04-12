from PyQt5 import QtWidgets, QtCore, QtGui
import copy
from grale.constants import * # TODO
from base import Layer, EmptyGraphicsItem, PointGraphicsItemBase

class SinglePointGraphicsItem(PointGraphicsItemBase):

    s_fontBrush = QtGui.QBrush(QtCore.Qt.green)

    #s_normalBrush = QtGui.QBrush(QtCore.Qt.yellow)
    s_normalBrush = QtGui.QBrush(QtGui.QColor(255, 255, 0, 128))
    #s_selectedBrush = QtGui.QBrush(QtCore.Qt.blue)
    s_selectedBrush = QtGui.QBrush(QtGui.QColor(0, 0, 255, 128))

    s_normalPen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.gray), 0.2)
    s_selectedPen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.darkGray), 0.2)

    def __init__(self, parent, uuid, timedelay, label):
        super(SinglePointGraphicsItem, self).__init__(parent, uuid, label)

        self.setFlag(QtWidgets.QGraphicsItem.ItemSendsScenePositionChanges)
        self.parent = parent

        self.circle = QtWidgets.QGraphicsEllipseItem(-1, -1, 2, 2, self)
        self.circle.setPen(SinglePointGraphicsItem.s_normalPen)
        self.circle.setBrush(SinglePointGraphicsItem.s_normalBrush)

        self.line1 = QtWidgets.QGraphicsLineItem(-0.2,-0.2, 0.2, 0.2, self)
        self.line1.setPen(SinglePointGraphicsItem.s_normalPen)

        self.line2 = QtWidgets.QGraphicsLineItem(-0.2, 0.2, 0.2, -0.2, self)
        self.line2.setPen(SinglePointGraphicsItem.s_normalPen)

        self.tdtxt = QtWidgets.QGraphicsSimpleTextItem(self)
        self.tdtxt.setPen(QtGui.QPen(SinglePointGraphicsItem.s_fontBrush, 0.2))
        self.tdtxt.setBrush(SinglePointGraphicsItem.s_fontBrush)

        if timedelay is not None:
            self._setTDText(timedelay)
            self.tdtxt.setVisible(True)
        else:
            self.tdtxt.setVisible(False)

    def _setTDText(self, timedelay):
        self.tdtxt.setText("{} days".format(timedelay))
        r = self.tdtxt.boundingRect()
        dx = (r.right()+r.left())/2.0
        self.tdtxt.setTransform(QtGui.QTransform(0.1,0,0,-0.1,0,0))
        self.tdtxt.setTransform(QtGui.QTransform(1,0,0,1,-dx, -20), True)

    def onSelected(self, value):
        if value:
            self.circle.setPen(SinglePointGraphicsItem.s_selectedPen)
            self.circle.setBrush(SinglePointGraphicsItem.s_selectedBrush)
        else:
            self.circle.setPen(SinglePointGraphicsItem.s_normalPen)
            self.circle.setBrush(SinglePointGraphicsItem.s_normalBrush)

    def itemChange(self, change, value):
        if change == QtWidgets.QGraphicsItem.ItemPositionHasChanged:
            self.parent.onPointPositionChanged(self.getUuid(), (value.x(), value.y()))

        return super(SinglePointGraphicsItem, self).itemChange(change, value)

    def fetchSettings(self):
        ptInfo = super(SinglePointGraphicsItem, self).fetchSettings()

        if ptInfo["timedelay"]:
            self.tdtxt.setVisible(True)
            self._setTDText(ptInfo["timedelay"])
        else:
            self.tdtxt.setVisible(False)

class TriangleItem(QtWidgets.QGraphicsPolygonItem):

    s_normalBrush = QtGui.QBrush(QtGui.QColor(192, 192, 192, 128))
    s_selectedBrush = QtGui.QBrush(QtGui.QColor(192, 192, 255, 192))

    s_normalPen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.darkGray), 2)
    s_selectedPen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.darkGray), 2)

    s_normalPen.setCosmetic(True)
    s_selectedPen.setCosmetic(True)

    def __init__(self, parent, uuid, ptsCoords):
        super(TriangleItem, self).__init__(parent)
        
        self.layer = parent.getLayer()
        self.uuid = uuid
        self.ptsCoords = ptsCoords
        self._updatePolygon()
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable)
        self.onSelected(False)

    def getLayer(self):
        return self.layer

    def getUuid(self):
        return self.uuid

    def updatePointPosition(self, ptKey, pos):
        self.ptsCoords[ptKey] = pos
        self._updatePolygon()

    def _updatePolygon(self):
        self.setPolygon(QtGui.QPolygonF([ QtCore.QPointF(self.ptsCoords[k][0], self.ptsCoords[k][1]) for k in self.ptsCoords ]))

    def paint(self, painter, option, widget):
        painter.setPen(self.pen())
        painter.setBrush(self.brush())
        painter.drawPolygon(self.polygon())

    def onSelected(self, value):
        if value:
            self.setBrush(TriangleItem.s_selectedBrush)
            self.setPen(TriangleItem.s_selectedPen)
        else:
            self.setBrush(TriangleItem.s_normalBrush)
            self.setPen(TriangleItem.s_normalPen)

    def itemChange(self, change, value):
        if change == QtWidgets.QGraphicsItem.ItemSelectedChange:
            self.onSelected(value)

        return super(TriangleItem, self).itemChange(change, value)

class MultiplePointGraphicsItem(EmptyGraphicsItem):
    def __init__(self, layer, parent = None):
        super(MultiplePointGraphicsItem, self).__init__(parent)

        self.layer = layer
        self.pointItems = { }
        self.triangleItems = { }

        self.updatePoints()

    def getLayer(self):
        return self.layer

    def addPoint(self, pos, timedelay = None, label = None, uuid = None):
        if not uuid:
            uuid = self.layer.addPoint(pos, timedelay=timedelay, label=label)
        else:
            self.layer.setPoint(uuid, pos, timedelay=timedelay, label=label)
        self.updatePoints()
        return uuid

    def _addPointItem(self, uuid, pos, timedelay, label):
        # TODO: check if uuid already present as safety check?

        m = SinglePointGraphicsItem(self, uuid, timedelay, label)
        m.setZValue(1)
        m.setPos(pos[0], pos[1])
        self.pointItems[uuid] = m

    def _addTriangItem(self, uuid, triangInfo):
        t = TriangleItem(self, uuid, triangInfo)
        t.setZValue(0)
        self.triangleItems[uuid] = t

    def updatePoints(self):
        
        # clear existing points
        s = self.scene()

        for itemDict in [ self.triangleItems, self.pointItems ]:
            for i in itemDict:
                if s: s.removeItem(itemDict[i])

        self.pointItems = { }
        self.triangleItems = { }

        pts = self.layer.getPoints()
        for k in pts:
            pos = pts[k]["xy"]
            self._addPointItem(k, pos, pts[k]["timedelay"], pts[k]["label"])

        triangs = self.layer.getTriangles()
        for k in triangs:
            t = { }
            for ptKey in triangs[k]:
                t[ptKey] = self.layer.getPoint(ptKey)["xy"]

            self._addTriangItem(k, t)
    
    def onPointPositionChanged(self, ptKey, newPos):
        triangs = self.layer.getTrianglesContainingPoint(ptKey)
        for tKey in triangs:
            if tKey in self.triangleItems: # Can be the case during initialization
                self.triangleItems[tKey].updatePointPosition(ptKey, newPos)


class PointsLayer(Layer):
    def __init__(self, name = "Points layer"):
        super(PointsLayer, self).__init__(name)

        self.triangles = { }
        self.trianglesInPoint = { } # For easier deletion

    def getTrianglesContainingPoint(self, ptKey):
        if not ptKey in self.trianglesInPoint:
            return []
        return copy.deepcopy(self.trianglesInPoint[ptKey])

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

            self.setPoint(key, (xy[0], xy[1]), timedelay=td, label=label)
            pointsAdded.add(key)

        if triangles and imgDat.hasTriangulation():
            # Keep track of the triangles a point's involved in. This helps
            # when deleting a point
            for triangKey in triangList:
                self.triangles[triangKey] = triangList[triangKey]
                trianglesAdded.add(triangKey)

                for ptKey in triangList[triangKey]:
                    if not ptKey in self.trianglesInPoint:
                        self.trianglesInPoint[ptKey] = set()

                    self.trianglesInPoint[ptKey].add(triangKey)
        
        return pointsAdded, trianglesAdded

    def _setPointKw(self, pointDict, timedelay, label):
        pointDict["timedelay"] = timedelay
        pointDict["label"] = label

    def clearPoint(self, key):
        #import pprint
        #print("Points:")
        #pprint.pprint(self.points)
        #print("Triangles:")
        #pprint.pprint(self.triangles)
        #print("Triangles in points:")
        #pprint.pprint(self.trianglesInPoint)

        super(PointsLayer, self).clearPoint(key)
        affectedTriangles = { }
        if key in self.trianglesInPoint:
            for triangKey in self.trianglesInPoint[key]:
                otherPoints = self.triangles[triangKey]
                for k in otherPoints:
                    if k != key and k in self.trianglesInPoint:
                        self.trianglesInPoint[k].remove(triangKey)

                affectedTriangles[triangKey] = self.triangles[triangKey]
                del self.triangles[triangKey]

            del self.trianglesInPoint[key]

        return affectedTriangles

    def clearAllPoints(self):
        super(PointsLayer, self).clearAllPoints()
        self.triangles = { }
        self.trianglesInPoint = { }

    def addTriangle(self, p0Key, p1Key, p2Key, triangKey = None):

        pts = set([ p0Key, p1Key, p2Key ])
        if len(pts) != 3:
            raise Exception("You need to specify three distinct points")

        # Check that the points exist
        for key in pts:
            self.getPoint(key)

        if triangKey:
            if triangKey in self.triangles:
                raise Exception("Specified triangle key already exists")
        else:
            triangKey = self.generateKey()
        self.triangles[triangKey] = pts

        for ptKey in pts:
            if not ptKey in self.trianglesInPoint:
                self.trianglesInPoint[ptKey] = set()
            self.trianglesInPoint[ptKey].add(triangKey)

        return triangKey

    def getTriangles(self):
        return copy.deepcopy(self.triangles)

    def getTriangle(self, key):
        return copy.deepcopy(self.triangles[key])

    def clearTriangle(self, key):
        pts = self.triangles[key]
        for p in pts:
            self.trianglesInPoint[p].remove(key)
        del self.triangles[key]

    def createGraphicsItem(self):
        return MultiplePointGraphicsItem(self)

    def toSettings(self):
        pts = self.getPoints()
        ptIdx = 0
        ptKeyIdx = { }
        for k in pts:
            ptKeyIdx[k] = ptIdx
            ptIdx += 1

        settings = {
                "name": self.getName(),
                "points": [ pts[k] for k in pts ],
                "triangles": [ list(map(lambda x: ptKeyIdx[x], self.triangles[t])) for t in self.triangles ]
        }

        return settings

    @staticmethod
    def fromSettings(settings):
        l = PointsLayer(settings["name"])
        keys = [ ]
        for p in settings["points"]:
            p = copy.deepcopy(p)
            xy = p["xy"]
            del p["xy"]
            keys.append(l.addPoint(xy, **p))

        for pts in settings["triangles"]:
            pts = list(map(lambda x: keys[x], pts))
            l.addTriangle(*pts)
            
        return l
