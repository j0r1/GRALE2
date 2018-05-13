from PyQt5 import QtWidgets, QtCore, QtGui
import numpy as np
import sys
import pprint
import uuid
import json
import copy

import grale.images as images # TODO
import grale.contourfinder as contourfinder

from base import GraphicsView, GraphicsScene, PointGraphicsItemBase, Layer, TriangleItem
from imagelayer import FITSImageLayer, RGBImageLayer, RGBGraphicsItem, FITSGraphicsItem, ImageGraphicsItem
from pointslayer import PointsLayer, MultiplePointGraphicsItem
from checkqt import checkQtAvailable
from hullprocessor import createTriangulationFromHull
from actionstack import ActionStack, objectToDouble
from contourleveldialog import ContourLevelDialog
from pointinfodialog import PointInfoDialog
from fitslayerinfodialog import FITSLayerInfoDialog

class DrawItem(QtWidgets.QGraphicsPathItem):
    def __init__(self, parent = None):
        super(DrawItem, self).__init__(parent)

        pen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.green), 2)
        pen.setCosmetic(True)

        self.setPen(pen)

class LayerScene(GraphicsScene):
    def __init__(self):
        super(LayerScene, self).__init__()

        self.actionStack = self.allocateActionStackInstance()

        self.drawItem = DrawItem()
        self.drawItem.setVisible(False)
        self.drawItem.setZValue(10000)
        self.addItem(self.drawItem)

        self.drawPath = [ ]

    def allocateActionStackInstance(self):
        return ActionStack()

    def getLayerItem(self, uuid):
        raise Exception("Implement in derived class!")

    def undo(self):
        self.actionStack.undo(self)

    def redo(self):
        self.actionStack.redo(self)

    def getActionStack(self):
        return self.actionStack

    def getNewPointLabel(self):
        return None

    def _addPoint(self, pos, label = None):
        curItem, curLayer = self.getCurrentItemAndLayer()
        if not curItem or not curLayer:
            return

        uuid = curItem.addPoint(pos[0], pos[1], label)
        ptInfo = curLayer.getPoint(uuid)
        if type(curLayer) == PointsLayer:
            self.getActionStack().recordAddNormalPoint(curLayer.getUuid(), uuid, ptInfo["xy"], ptInfo["timedelay"], ptInfo["label"], None)
        else:
            self.getActionStack().recordAddMatchPoint(curLayer.getUuid(), uuid, ptInfo["xy"], ptInfo["label"], None)

    def onPointLabelChanged(self, layerUuid, pointUuid, oldLabel, newLabel):
        self.getActionStack().recordLabelChange(layerUuid, pointUuid, oldLabel, newLabel)

    def mouseHandler_moved(self, movableItem, cursorPos, itemPos, moveId):
        if movableItem:
            oldPosition, newPosition = movableItem.syncPosition()
            layerUuid = movableItem.getLayer().getUuid()
            self.getActionStack().recordPointMove(layerUuid, movableItem.getUuid(), oldPosition, newPosition, moveId)

    def mouseHandler_click(self, movableItem, pos, modInfo):
        selItems = self.selectedItems() # Note, the following 'super' might clear the selection
        super(LayerScene, self).mouseHandler_click(movableItem, pos, modInfo)

        if movableItem is None: # possibly new point
            # If nothing is selected yet, add a point
            if len(selItems) == 0: 
                fi = self.focusItem()
                if not fi:
                    label = self.getNewPointLabel()
                    self._addPoint(pos, label)
                else:
                    fi.clearFocus()
            else: # clear selection and focus
                for i in selItems:
                    i.setSelected(False)

                fi = self.focusItem()
                if fi:
                    fi.clearFocus()

    def _checkUniqueLabels(self, matchPoints, identifier):
        labels = set()
        for uuid in matchPoints:
            l = matchPoints[uuid]["label"]
            if l in labels:
                raise Exception("Label '{}' is present multiple times in {}".format(l, identifier))

            labels.add(l)

    def _matchTransform(self, curLayer, curItem, refLayer):
        rgbPoints = curLayer.getPoints(True)  # True = transformed to scene coords
        fitsPoints = refLayer.getPoints(True)
        try:
            self._checkUniqueLabels(rgbPoints, "RGB image")
            self._checkUniqueLabels(fitsPoints, "FITS image")
        except Exception as e:
            self.warning("Error checking unique labels", str(e))
            return

        matchPoints = { }

        # Find points to match
        for pts in [ rgbPoints, fitsPoints ]:
            for uuid in pts:
                label = pts[uuid]["label"]
                xy = pts[uuid]["xy"]

                if not label in matchPoints:
                    matchPoints[label] = [ xy ]
                else:
                    matchPoints[label].append(xy)

        # only use matching if there are two entries, first is rgb, second is fits
        pointsToMatch = [ matchPoints[label] for label in matchPoints if len(matchPoints[label]) == 2 ]
        
        try:
            oldTrans = curLayer.getImageTransform()
            curLayer.matchToPoints(pointsToMatch)
            newTrans = curLayer.getImageTransform()
        except Exception as e:
            self.warning("Can't calculate correspondence", str(e))
            return

        self.getActionStack().recordRGBTransform(curLayer.getUuid(), oldTrans, newTrans)
        curItem.updateTransform()

    def _getInteractiveContour(self, pos):
        # store original visibilities
        itemList = [ i for i in self.items() if type(i) == MultiplePointGraphicsItem ] + [ self.drawItem, self.getAxesItem() ]
        itemList = [ (i,i.isVisible()) for i in itemList ]
        for i,v in itemList:
            i.setVisible(False)

        try:
            pixelSize = 256 # TODO: make this configurable
            viewSize = 4 # 4 arcsec, TODO: make this configurable
            pixelBlurSize = 4 # TODO: make this configurable

            bottomLeft, topRight = [pos[0] - viewSize/2.0, pos[1] - viewSize/2.0], [pos[0] + viewSize/2.0, pos[1] + viewSize/2.0]
            arr = self.getSceneRegionNumPyArray(bottomLeft, topRight, pixelSize, pixelSize, grayScale=True)
            arr = arr.astype(np.double)
            
            # Blur the region to get smoother contours
            from scipy.ndimage.filters import gaussian_filter
            arr = gaussian_filter(arr, sigma=pixelBlurSize)
        finally:
            # Restore original visibilities
            for i,v in itemList:
                i.setVisible(v)
        
        # From https://stackoverflow.com/a/21613346/2828217
        def isInside(pos, cntr):
            from shapely.geometry import Point
            from shapely.geometry.polygon import Polygon

            try:
                point = Point(pos[0], pos[1])
                polygon = Polygon(cntr)
                return polygon.contains(point)
            except Exception as e:
                print("WARNING: {}".format(e))
                return False

        levels = list(np.linspace(arr.min(),arr.max(), 255))
        cf = contourfinder.ContourFinder(arr, bottomLeft, topRight)
        contours = cf.findMultipleContours(levels)
        filteredContours = [ ]
        for c in contours:
            for part in c:
                if isInside(pos, part):
                    filteredContours.append(part)
                    break
    
        contours = filteredContours
        if len(contours) < 1:
            return None

        self.drawItem.setVisible(True)
        dlg = ContourLevelDialog(contours, self.drawItem, self.getDialogWidget())
        r = contours[dlg.getSelectedContour()] if dlg.exec_() else None
        self.drawItem.setVisible(False)

        return r

    def _recenter(self, item, layer, pos):
        pixPos = item.getPixelPosition(pos)
        ra, dec = layer.getPixelCoordinates(pixPos[0], pixPos[1])
        oldRa, oldDec = layer.getCenter()
        layer.setCenter(ra, dec)
        item.updateTransform()
        self.getActionStack().recordCenterLocation(layer.getUuid(), [oldRa, oldDec], [ra, dec])

    def mouseHandler_doubleClick(self, movableItem, pos, modInfo):
        if not movableItem:
            item, layer = self.getCurrentItemAndLayer()
            if type(layer) == FITSImageLayer:
                self._recenter(item, layer, pos)
            elif type(layer) == RGBImageLayer: # Point match
                refLayer = self.getPointMatchReferenceLayer()
                if refLayer:
                    self._matchTransform(layer, item, refLayer)
            elif type(layer) == PointsLayer:
                contour = self._getInteractiveContour(pos)
                if contour is not None:
                    self._addTriangulationFromHull(contour)
        else:
            if not movableItem.isNormalPoint():
                movableItem.toggleFocus()
            else: # Point item
                self._pointDoubleClicked(movableItem)

    def _pointDoubleClicked(self, movableItem):
        if movableItem.isSelected():
            pts = [ ]
            for i in self.selectedItems():
                i = PointGraphicsItemBase.getPointGraphicsItem(i) 
                if i and i.isNormalPoint():
                    pts.append(i)

            if len(pts) < 3:
                movableItem.toggleFocus()
                return
            
            layer = pts[0].getLayer()
            for i in pts:
                if i.getLayer() is not layer:
                    self.warning("Can't create triangulation", "Not all selected points are in the same layer")
                    return

            layerItem = self.getLayerItem(layer.getUuid())

            ids = [ i.getUuid() for i in pts ]
            coords = [ [ i.pos().x(), i.pos().y() ] for i in pts ]

            from scipy.spatial import Delaunay

            try:
                if len(pts) < 3:
                    raise Exception("Need at least three points")

                tri = Delaunay(coords)
            except Exception as e:
                self.warning("Can't create triangulation", "Can't create triangulation: {}".format(e))
                return
            
            addId = uuid.uuid4()
            layerUuid = layer.getUuid()

            for pts in tri.simplices:
                ptIds = [ ids[p] for p in pts ]
                tUuid = layerItem.addTriangle(*ptIds)
                self.getActionStack().recordAddTriangle(layerUuid, tUuid, ptIds, addId)

        else:
            movableItem.toggleFocus()

    def _deleteSelectedItems(self, deletePoints, deleteTriangles):
        items = self.selectedItems()
        pointItems, triangItems = [ ], [ ]
        for i in items:
            p = PointGraphicsItemBase.getPointGraphicsItem(i)
            if p:
                pointItems.append(p)
            else:
                t = TriangleItem.getTriangleItem(i)
                if t:
                    triangItems.append(t)

        self._deleteItems(pointItems, triangItems, deletePoints, deleteTriangles)

    def _deleteItems(self, pointItems, triangItems, deletePointsFlag, deleteTrianglesFlag):
    
        affectedLayerItems = set()
        deleteId = uuid.uuid4()

        deleteNormalPoints = [ ]
        deleteTriangles = [ ]
        deleteMatchPoints = [ ]

        needUpdate = False
        if triangItems and deleteTrianglesFlag:
            for t in triangItems:
                tUuid = t.getUuid()
                layer = t.getLayer()
                layerUuid = layer.getUuid() # TODO: better way to get layer item from layer!
                layerItem = self.getLayerItem(layer.getUuid())
                pts = layer.getTriangle(tUuid)

                #self.getActionStack().recordDeleteTriangle(layerUuid, tUuid, pts, deleteId)
                deleteTriangles.append({"layer": layerUuid, "triangle": tUuid, "points": pts })

                layerItem.clearTriangle(tUuid)
        
        if pointItems and deletePointsFlag:
            for p in pointItems:
                pUuid = p.getUuid()
                layer = p.getLayer()
                layerUuid = layer.getUuid()
                layerItem = self.getLayerItem(layer.getUuid())
                ptInfo = layer.getPoint(pUuid)

                isNormalPoint = p.isNormalPoint() # Note: we need to do this before deleting the point, as this will also delete the graphics item
                triangs = layerItem.clearPoint(pUuid) # returns affected triangles

                if isNormalPoint:
                    deleteNormalPoints.append({"layer": layerUuid, "point": pUuid, "pos": ptInfo["xy"], "timedelay": ptInfo["timedelay"], "label": ptInfo["label"] })
                    for t in triangs:
                        deleteTriangles.append({"layer": layerUuid, "triangle": t, "points": triangs[t] })
                else: # Match point
                    pos = ptInfo["xy"]
                    pos = layer.getImageTransform().map(pos[0], pos[1])
                    deleteMatchPoints.append({"layer": layerUuid, "point": pUuid, "pos": pos, "label": ptInfo["label"] })

        if deleteNormalPoints:
            self.getActionStack().recordDeletePoints(deleteNormalPoints, deleteId)
        if deleteTriangles:
            self.getActionStack().recordDeleteTriangles(deleteTriangles, deleteId)
        if deleteMatchPoints:
            self.getActionStack().recordDeleteMatchPoints(deleteMatchPoints, deleteId)

    def _addTriangulationFromHull(self, path):

        try:
            pts, simp = createTriangulationFromHull(np.array(path)) # TODO: default settings
        except Exception as e:
            self.warning("Can't create triangulation", "Can't create hull points and triangulation: {}".format(e))
            return

        addId = uuid.uuid4()
        item, layer = self.getCurrentItemAndLayer()
        layerUuid = layer.getUuid()

        if type(layer) != PointsLayer:
            self.warning("Can't create triangulation", "Current layer is not a points layer, can't create a triangulation based on the drawn hull")
            return

        uuids = [ ]

        importPoints = [ ]
        for p in pts:
            ptUuid = item.addPoint(float(path[p][0]), float(path[p][1]))
            uuids.append(ptUuid)

            ptInfo = layer.getPoint(ptUuid)
            importPoints.append({ "layer": layerUuid, "point": ptUuid, "pos": ptInfo["xy"], "timedelay": ptInfo["timedelay"], "label": ptInfo["label"] })
        
        self.getActionStack().recordAddNormalPoints(importPoints, addId)

        importTriangles = [ ]
        for p0, p1, p2 in simp:
            tUuid = item.addTriangle(uuids[p0], uuids[p1], uuids[p2])
            importTriangles.append({"triangle": tUuid, "layer": layerUuid, "points": [ uuids[p0], uuids[p1], uuids[p2] ] })
        
        self.getActionStack().recordAddTriangles(importTriangles, addId)
        

    def keyHandler_clicked(self, key, keyStr, modifiers):

        if keyStr == "Z" and modifiers["control"] == True and modifiers["shift"] == False:
            self.undo()
            return

        if keyStr == "Z" and modifiers["control"] == True and modifiers["shift"] == True:
            self.redo()
            return

        if keyStr == "C" and modifiers["control"] == True and modifiers["shift"] == False and modifiers["alt"] == False:
            self.copy()
            return

        if keyStr == "V" and modifiers["control"] == True and modifiers["shift"] == False and modifiers["alt"] == False:
            self.paste()
            return

        if keyStr == "X" and modifiers["control"] == True and modifiers["shift"] == False and modifiers["alt"] == False:
            self.cut()
            return

        if key == QtCore.Qt.Key_Delete or key == QtCore.Qt.Key_Backspace:
            self._deleteSelectedItems(modifiers["shift"], modifiers["control"])
            return

        if keyStr == "I":
            for i in self.items():
                print(i)

    def importFromImagesData(self, imgDat, imgIdx):
        item, layer = self.getCurrentItemAndLayer()
        if type(layer) != PointsLayer:
            self.warning("Can't import images data", "Can't import images data: the current layer is not a points layer")
            return

        layerUuid = layer.getUuid()
        pts, triangs = layer.importFromImagesData(imgDat, imgIdx) # TODO: specify flags?
        addId = uuid.uuid4()

        importPoints = [ ]
        for ptUuid in pts:
            ptInfo = layer.getPoint(ptUuid)
            importPoints.append({ "layer": layerUuid, "point": ptUuid, "pos": ptInfo["xy"], "timedelay": ptInfo["timedelay"], "label": ptInfo["label"] })
        
        self.getActionStack().recordAddNormalPoints(importPoints, addId)

        importTriangles = [ ]
        for tUuid in triangs:
            uuids = layer.getTriangle(tUuid)
            importTriangles.append({"triangle": tUuid, "layer": layerUuid, "points": uuids })
        
        self.getActionStack().recordAddTriangles(importTriangles, addId)

        item.updatePointList([u for u in pts], [t for t in triangs])

    def mouseHandler_startDrawing(self, pos, modInfo):
        item, layer = self.getCurrentItemAndLayer()
        if type(layer) != PointsLayer:
            return False

        #pprint.pprint(modInfo) 

        self.drawItem.setVisible(True)
        path = QtGui.QPainterPath(QtCore.QPointF(pos[0], pos[1]))
        self.drawItem.setPath(path)

        self.drawPath = [ pos ]
        return True

    def mouseHandler_drawing(self, pos):
        path = self.drawItem.path()
        path.lineTo(pos[0], pos[1])
        self.drawItem.setPath(path)

        self.drawPath.append(pos)

    def mouseHandler_stopDrawing(self, pos):
        self.drawItem.setVisible(False)

        self.drawPath.append(pos)
        path = self.drawPath
        self.drawPath = [ ]

        if len(path) < 10: # probably not intended, just ignore
            return
 
        self._addTriangulationFromHull(path)

    def copy(self):
        items = self.selectedItems()
        pointItems, triangItems, matchPointItems = [ ], [ ], [ ]
        for i in items:
            p = PointGraphicsItemBase.getPointGraphicsItem(i)
            if p:
                if p.isNormalPoint():
                    pointItems.append(p)
                else:
                    matchPointItems.append(p)
            else:
                t = TriangleItem.getTriangleItem(i)
                if t:
                    triangItems.append(t)

        ptMap = { }

        copyInfo = {
                "normalpoints": [ None ] * len(pointItems),
                "matchpoints": [ None ] * len(matchPointItems),
                "triangles": [ None ] * len(triangItems)
        }

        for k in range(len(pointItems)):
            i = pointItems[k]
            uuid, layer = i.getUuid(), i.getLayer()
            ptMap[uuid] = k
            copyInfo["normalpoints"][k] = layer.getPoint(uuid)

        # It's possible that some triangle corners were not selected,
        # take these into account as well
        for i in triangItems:
            layer = i.getLayer()
            pts = layer.getTriangle(i.getUuid())
            for uuid in pts:
                if not uuid in ptMap:
                    k = len(copyInfo["normalpoints"])
                    copyInfo["normalpoints"].append(layer.getPoint(uuid))
                    ptMap[uuid] = k

        for k in range(len(triangItems)):
            i = triangItems[k]
            uuid, layer = i.getUuid(), i.getLayer()
            ptList = list(map(lambda x: ptMap[x], layer.getTriangle(uuid)))
            copyInfo["triangles"][k] = ptList

        for k in range(len(matchPointItems)):
            i = matchPointItems[k]
            uuid, layer = i.getUuid(), i.getLayer()
            layerItem = self.getLayerItem(layer.getUuid())
            ptInfo = layer.getPoint(uuid)
            ptInfo["xy"] = layerItem.getScenePosition(ptInfo["xy"]) # make sure the scene coords are copied
            copyInfo["matchpoints"][k] = ptInfo

        copyStr = json.dumps(copyInfo, sort_keys=True, indent=4, separators=(',', ': '))
        clipboard = QtWidgets.QApplication.clipboard()
        clipboard.setText(copyStr)

        return pointItems + matchPointItems, triangItems

    def paste(self):
        clipboard = QtWidgets.QApplication.clipboard()
        txt = clipboard.text()
        try:
            copyInfo = json.loads(txt)
            item, layer = self.getCurrentItemAndLayer()
            layerUuid = layer.getUuid()
    
            addedNormalPoints = [ ]
            addedTriangles = [ ]
            addedMatchPoints = [ ]

            addId = uuid.uuid4()
            if type(layer) == PointsLayer:
                ptMap = { }
                for pIdx in range(len(copyInfo["normalpoints"])):
                    added = True
                    
                    ptInfo = copyInfo["normalpoints"][pIdx]
                    pos = ptInfo["xy"]

                    ptUuid = item.addPoint(pos[0], pos[1], ptInfo["label"], objectToDouble(ptInfo["timedelay"]))
                    ptMap[pIdx] = ptUuid
                    
                    addedNormalPoints.append({"layer": layerUuid, "point": ptUuid, "pos": pos, "timedelay": ptInfo["timedelay"], "label": ptInfo["label"] })

                for tIdx in range(len(copyInfo["triangles"])):
                    pts = list(map(lambda x: ptMap[x], copyInfo["triangles"][tIdx]))
                    tUuid = item.addTriangle(*pts)
                    addedTriangles.append({"layer": layerUuid, "triangle": tUuid, "points": pts })

            else: # image layer
                for pIdx in range(len(copyInfo["matchpoints"])):
                    ptInfo = copyInfo["matchpoints"][pIdx]

                    pos = ptInfo["xy"]
                    ptUuid = item.addPoint(pos[0], pos[1], ptInfo["label"])
                    addedMatchPoints.append({"layer": layerUuid, "point": ptUuid, "pos": pos, "label": ptInfo["label"]})
            
            if addedNormalPoints:
                self.getActionStack().recordAddNormalPoints(addedNormalPoints, addId)
            if addedTriangles:
                self.getActionStack().recordAddTriangles(addedTriangles, addId)
            if addedMatchPoints:
                self.getActionStack().recordAddMatchPoints(addedMatchPoints, addId)

        except Exception as e:
            self.warning("Unable to paste clipboard contents", "Can't paste: {}".format(e))

    def cut(self):
        ptItems, triangItems = self.copy()
        self._deleteItems(ptItems, triangItems, True, True)

    def _editPointInfo(self, pointItem):
        layer = pointItem.getLayer()
        ptUuid = pointItem.getUuid()

        ptInfo = layer.getPoint(ptUuid)
        dlg = PointInfoDialog(ptInfo, self.getDialogWidget())
        if dlg.exec_():
            ptInfoNew = dlg.getPointInfo()
            layer.setPoint(ptUuid, ptInfoNew["xy"], ptInfoNew["label"], ptInfoNew["timedelay"])
            pointItem.fetchSettings() # adjust the QGraphicsItem to match the point settings

            self.getActionStack().recordSetNormalPoint(layer.getUuid(), ptUuid, ptInfo["xy"], ptInfoNew["xy"],
                                      ptInfo["timedelay"], ptInfoNew["timedelay"], ptInfo["label"], ptInfoNew["label"])


    def _editFITSLayerInfo(self, fitsItem):
        layer = fitsItem.getLayer()
        oldMinMax = tuple(map(float, layer.getMinMax()))
        oldRaDec = tuple(map(float, layer.getCenter()))
        dlg = FITSLayerInfoDialog(oldRaDec, oldMinMax, self.getDialogWidget())
        if dlg.exec_():
            newMinMax = dlg.getMinMax()
            newRaDec = dlg.getCenter()

            if newMinMax == oldMinMax: # only a change in center
                layer.setCenter(*newRaDec)
                fitsItem.updateTransform()
                self.getActionStack().recordCenterLocation(layer.getUuid(), oldRaDec, newRaDec)
            else:
                layer.setCenter(*newRaDec)
                layer.setMinMax(*newMinMax)
                fitsItem.updatePixmap()
                fitsItem.updateTransform()
                self.getActionStack().recordCenterAndMinMax(layer.getUuid(), oldRaDec, newRaDec, oldMinMax, newMinMax)

    def mouseHandler_rightClick(self, pointItem, pos, modInfo):
        if pointItem:
            pointItem = PointGraphicsItemBase.getPointGraphicsItem(pointItem)
            if pointItem and pointItem.isNormalPoint():
                self._editPointInfo(pointItem)
        else:
            # Look for RGBGraphicsItem or FITSGraphicsItem
            for item in self.items(QtCore.QPointF(pos[0], pos[1])):
                if type(item) == FITSGraphicsItem:
                    self._editFITSLayerInfo(item)
                    break


    def getCurrentItemAndLayer(self):
        raise Exception("Implement in derived class")

    def getPointMatchReferenceLayer(self):
        pass

class SingleLayerScene(LayerScene):
    def __init__(self, layer):
        super(SingleLayerScene, self).__init__()

        self.layer = layer
        self.item = layer.createGraphicsItem()
        self.addItem(self.item)

    def getLayerItem(self, uuid):
        if uuid == self.layer.getUuid():
            return self.item
        raise Exception("Layer with uuid {} not found".format(uuid))

    def getCurrentItemAndLayer(self):
        return self.item, self.layer

class PointsSingleLayerScene(SingleLayerScene):
    def __init__(self, layer, bgLayer = None):
        super(PointsSingleLayerScene, self).__init__(layer)

        curItem, curLayer = self.getCurrentItemAndLayer()
        curItem.setZValue(1)

        if bgLayer:
            self.bgItem = bgLayer.createGraphicsItem()
            self.bgItem.setZValue(0)
            self.bgItem.showPoints(False)
            self.addItem(self.bgItem)

class ImageSingleLayerScene(SingleLayerScene):
    def __init__(self, layer):
        super(ImageSingleLayerScene, self).__init__(layer)

        self.counter = 0

    def getNewPointLabel(self):
        self.counter += 1
        return str(self.counter)

class FITSSingleLayerScene(ImageSingleLayerScene):
    def __init__(self, layer):
        super(FITSSingleLayerScene, self).__init__(layer)

class MatchRGBLayerScene(ImageSingleLayerScene):
    def __init__(self, rgbLayer, matchLayer):
        super(MatchRGBLayerScene, self).__init__(rgbLayer)

        self.matchLayer = matchLayer

    def getPointMatchReferenceLayer(self):
        return self.matchLayer


def matchTest():

    app = QtWidgets.QApplication(sys.argv)
    
    testFileName = "/home/jori/projects/graleeditor-hg/src/hst_12817_01_acs_wfc_f606w_drz.fits"
    fitsLayer = FITSImageLayer(testFileName)

    testFileName = "/tmp/040518_LG_dark-matter-galaxy_feat.jpg"
    #testFileName = "/home/jori/projects/a3827/eso1514a.jpg"
    rgbLayer = RGBImageLayer(testFileName)

    view1 = GraphicsView(FITSSingleLayerScene(fitsLayer))
    view2 = GraphicsView(MatchRGBLayerScene(rgbLayer, fitsLayer))

    status = app.exec_()
    view1.setScene(None) # To avoid a segfault
    view2.setScene(None) # To avoid a segfault

    pprint.pprint(fitsLayer.toSettings())
    pprint.pprint(rgbLayer.toSettings())

    return fitsLayer, rgbLayer

def pointsTest():
    import grale.images as images

    app = QtWidgets.QApplication(sys.argv)
    
    pointsLayer = PointsLayer()
    scene = PointsSingleLayerScene(pointsLayer)
    view = GraphicsView(scene)

    #scene.setAxisLeft(False)
    #scene.setPointSizeFixed(False)
    scene.importFromImagesData(images.ImagesData.load("/tmp/testimg.imgdata"), 0)
    #scene.importFromImagesData(images.ImagesData.load("/tmp/splittest.imgdat"), "all")

    status = app.exec_()
    
    view.setScene(None)

    pprint.pprint(pointsLayer.toSettings())
    return pointsLayer

def bgPointsTest():
    fitsLayer, rgbLayer = matchTest()
    
    app = QtWidgets.QApplication(sys.argv)

    pointsLayer = PointsLayer()
    view = GraphicsView(PointsSingleLayerScene(pointsLayer, rgbLayer))
    status = app.exec_()
    
    view.setScene(None)

    return pointsLayer

def matchSettingsTest():
    fitsLayer, rgbLayer = matchTest()

    app = QtWidgets.QApplication(sys.argv)

    newFits = Layer.fromSettings(fitsLayer.toSettings())
    newRgb = Layer.fromSettings(rgbLayer.toSettings())
    view1 = GraphicsView(FITSSingleLayerScene(newFits))
    view2 = GraphicsView(MatchRGBLayerScene(newRgb, newFits))
    status = app.exec_()
    
    view1.setScene(None)
    view2.setScene(None)

def main():
    checkQtAvailable()
    pointsTest()
    #matchTest()
    #bgPointsTest()
    
    #matchSettingsTest()

if __name__ == "__main__":
    main()
