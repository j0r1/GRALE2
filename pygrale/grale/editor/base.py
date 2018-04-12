from PyQt5 import QtWidgets, QtCore, QtGui
import copy
import uuid

# Layer can have points (normal points, match points), possibly with other
# properties (provided by derived class) (e.g. label for point matching,
# time delay information. The 'transform' is the transformation needed to
# convert the point coordinates to scene coordinates, which is mainly for
# point matching: there, the coordinates are in pixel space and the 
# transformation converts them to scene space (where each unit is an arc second)

class Layer(object):
    def __init__(self, name = "Untitled"):
        self.points = { }
        self.transform = QtGui.QTransform(1, 0, 0, 1, 0, 0)
        self.uuid = self.generateKey()
        self.name = name
    
    def setName(self, n):
        self.name = n

    def getName(self):
        return self.name

    def getUuid(self):
        return self.uuid
    
    def createGraphicsItem(self):
        raise Exception("Implement this in derived class!")

    def _setPointKw(self, **kwargs): # Implement this in derived class for other props
        pass

    def generateKey(self):
        key = uuid.uuid4()
        return key

    def addPoint(self, xy, **kwargs):
        key = self.generateKey()
        d = { "xy": [xy[0], xy[1]] }
        self._setPointKw(d, **kwargs)
        self.points[key] = d
        return key

    def setPoint(self, key, xy, **kwargs):
        d = { "xy": [xy[0], xy[1]] }
        self._setPointKw(d, **kwargs)
        self.points[key] = d

    def clearPoint(self, key):
        del self.points[key]

    def clearAllPoints(self):
        self.points = { }

    def getPoint(self, key):
        return copy.deepcopy(self.points[key])

    def getPoints(self, transformed = False):
        if not transformed:
            return copy.deepcopy(self.points)

        tranformedPoints = { }
        for key in self.points:
            newPt = copy.deepcopy(self.points[key]) # keep other properties
            xy = newPt["xy"]
            mXy = self.transform.map(QtCore.QPointF(xy[0], xy[1]))
            newPt["xy"] = [ mXy.x(), mXy.y() ]
            tranformedPoints[key] = newPt

        return tranformedPoints

    def getTransform(self):
        return self.transform

    def setTransform(self, t):
        self.transform = t

    @staticmethod
    def fromSettings(settings):
        from pointslayer import PointsLayer
        from imagelayer import RGBImageLayer, FITSImageLayer

        if "rgbfile" in settings:
            return RGBImageLayer.fromSettings(settings)
        if "fitsfile" in settings:
            return FITSImageLayer.fromSettings(settings)
        return PointsLayer.fromSettings(settings)

class EmptyGraphicsItem(QtWidgets.QGraphicsItem):
    def __init__(self, parent = None):
        super(EmptyGraphicsItem, self).__init__(parent)

    def boundingRect(self):
        # The childrenBoundingRect doesn't take the visibility status into account
        r = QtCore.QRectF(0,0,1,1)
        for i in self.childItems():
            if i.isVisible():
                r |= i.mapToParent(i.boundingRect()).boundingRect()
        return r

    def paint(self, painter, option, widget):
        pass

class SyncTextEditItem(QtWidgets.QGraphicsTextItem):
    def __init__(self, layer, uuid, label, parent):
        super(SyncTextEditItem, self).__init__(label, parent)

        self.layer = layer
        self.uuid = uuid

    def keyReleaseEvent(self, evt):
        #print("keyReleaseEvent")
        super(SyncTextEditItem, self).keyReleaseEvent(evt)
        self._syncValueAndCenter()
    
    def inputMethodEvent(self, evt):
        #print("inputMethodEvent")
        super(SyncTextEditItem, self).inputMethodEvent(evt)
        self._syncValueAndCenter()

    def focusOutEvent(self, evt):
        #print("focusOutEvent")
        super(SyncTextEditItem, self).focusOutEvent(evt)
        self._syncValueAndCenter()
        if not self.layer.getPoint(self.uuid)["label"]:
            self.setVisible(False)

    def recenter(self):
        r = self.boundingRect()
        dx = (r.right()+r.left())/2.0
        self.setTransform(QtGui.QTransform(0.1,0,0,-0.1,0,0))
        self.setTransform(QtGui.QTransform(1,0,0,1,-dx, 10), True)

    def _syncValueAndCenter(self):
        # Sync value
        label = self.document().toPlainText()
        ptInfo = self.layer.getPoint(self.uuid)
        xy, oldLabel = ptInfo["xy"], ptInfo["label"]
        del ptInfo["xy"]
        del ptInfo["label"]
        self.layer.setPoint(self.uuid, xy, label=label, **ptInfo)

        scene = self.scene()
        if scene:
            scene.onPointLabelChanged(self.layer.getUuid(), self.uuid, oldLabel, label)

        # Center
        self.recenter()

class LayerObjectGraphicsItem(EmptyGraphicsItem):
    def __init__(self, layer, objectUuid, parent):
        super(LayerObjectGraphicsItem, self).__init__(parent)
        self.layer = layer
        self.uuid = objectUuid

    def getUuid(self):
        return self.uuid

    def getLayer(self):
        return self.layer

class PointGraphicsItemBase(LayerObjectGraphicsItem):
    def __init__(self, parent, uuid, label = None):
        super(PointGraphicsItemBase, self).__init__(parent.getLayer(), uuid, parent)

        self.setMovable()

        label = "" if not label else label
        self.txt = SyncTextEditItem(self.getLayer(), uuid, label, self)
        self.txt.setFlag(QtWidgets.QGraphicsItem.ItemIsFocusable)
        self.txt.setDefaultTextColor(QtGui.QColor(QtCore.Qt.cyan))
        self.txt.setTextInteractionFlags(QtCore.Qt.TextEditorInteraction)
        self.txt.recenter()

        if not label:
            self.txt.setVisible(False)

        s = self.scene()
        if s:
            self.setTransform(s._getPointTransform())

    def isMovable(self):
        f = self.flags()
        if f&QtWidgets.QGraphicsItem.ItemIsMovable and f&QtWidgets.QGraphicsItem.ItemIsSelectable:
            return True
        return False

    def setMovable(self, v = True):
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, v)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, v)

    # Override in derived class to convert scene position to point coordinates
    def scenePositionToPointPosition(self, pos):
        return pos
    def pointPositionToScenePosition(self, pos):
        return pos

    def fetchSettings(self):
        layer = self.getLayer()
        props = layer.getPoint(self.uuid)
        pos = self.pointPositionToScenePosition(props["xy"])
       
        self.setPos(QtCore.QPointF(pos[0], pos[1]))

        if not props["label"]:
            self.txt.setVisible(False)
        else:
            self.txt.setPlainText(props["label"])
            self.txt.setVisible(True)
            self.txt.recenter()
        
        return props

    def syncPosition(self, pos = None):
        if not pos:
            pos = self.pos()
            pos = (pos.x(), pos.y())
        else:
            self.setPos(pos[0], pos[1])

        layer = self.getLayer()
        pos = self.scenePositionToPointPosition(pos)
        props = layer.getPoint(self.uuid)
       
        oldPosition, newPosition = props["xy"], pos
        del props["xy"]
        
        layer.setPoint(self.uuid, pos, **props)

        return oldPosition, newPosition

    def onSelected(self, value):
        print("TODO: implement onSelected in derived class")

    def itemChange(self, change, value):
        if change == QtWidgets.QGraphicsItem.ItemSelectedChange:
            self.onSelected(value)

        return super(PointGraphicsItemBase, self).itemChange(change, value)

    def toggleFocus(self):
        if self.txt.hasFocus():
            self.txt.clearFocus()
        else:
            self.txt.setFocus()
            self.txt.setVisible(True)

class GraphicsView(QtWidgets.QGraphicsView):
    def __init__(self, scene, initialZoom = 8.0, parent = None):
        super(GraphicsView, self).__init__(parent)

        self.setMouseTracking(True)

        self.setViewportUpdateMode(QtWidgets.QGraphicsView.BoundingRectViewportUpdate)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.setDragMode(QtWidgets.QGraphicsView.RubberBandDrag)

        self.setBackgroundBrush(QtCore.Qt.black)
        self.setScene(scene)

        self.setTransform(QtGui.QTransform(initialZoom, 0, 0, initialZoom, 0, 0))
        self.updateAxisOrientation()

        self.show()

    def updateAxisOrientation(self):
        zoom = self.getScale()
        scene = self.scene()
        x = scene.getAxisDirection()
        self.setTransform(QtGui.QTransform(x*zoom, 0, 0, -zoom, 0, 0), False)
        scene.onScaleChanged(zoom)

    def getCenter(self):
        return self.mapToScene(self.viewport().rect().center())

    def wheelEvent(self, event):
        if event.modifiers()&QtCore.Qt.ControlModifier:
            event.accept()

            #print(event.angleDelta())
            point = event.pos()
            pointPos = self.mapToScene(point)

            scaleFactor = 2**(-event.angleDelta().y() / 1200.0)
            self.scaleView(scaleFactor)

            newPointPos = self.mapFromScene(pointPos)
            diff = newPointPos-point
            diffScene = self.mapToScene(diff) - self.mapToScene(QtCore.QPoint(0,0))
            currentCenter = self.getCenter()

            self.centerOn(currentCenter+diffScene);
        else:
            super(GraphicsView, self).wheelEvent(event)

    def getScale(self):
        factor = self.transform().mapRect(QtCore.QRectF(0, 0, 1, 1)).width()
        return factor

    def setScale(self, factor):
        if factor < 0.01 or factor > 10000:
            return

        x = self.scene().getAxisDirection()
        t = QtGui.QTransform(x*factor, 0, 0, -factor, 0, 0)
        self.setTransform(t, False)

    def scaleView(self, scaleFactor):
        factor = self.transform().scale(scaleFactor, scaleFactor).mapRect(QtCore.QRectF(0, 0, 1, 1)).width()
        if factor < 0.01 or factor > 10000:
            return

        self.scale(scaleFactor, scaleFactor)

        self.scene().onScaleChanged(self.getScale())

        # TODO:emit signalZoomLevel(factor);

    def resizeEvent(self, evt):
        oldW, oldH = evt.oldSize().width(), evt.oldSize().height()
        if oldW < 1 or oldH < 1:
            return

        centerPos = self.mapToScene(QtCore.QPoint(oldW//2, oldH//2))
        self.centerOn(centerPos)


class AxesGraphicsItem(QtWidgets.QGraphicsItemGroup):
    def __init__(self, parent):
        super(AxesGraphicsItem, self).__init__(parent)
        linePen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.red), 2)
        linePen.setCosmetic(True)
        self.line1 = QtWidgets.QGraphicsLineItem(-1000,0, 1000, 0, self)
        self.line1.setPen(linePen)
        self.line2 = QtWidgets.QGraphicsLineItem(0,-1000, 0, 1000, self)
        self.line2.setPen(linePen)
        self.setZValue(1000000)

class GraphicsScene(QtWidgets.QGraphicsScene):

    def __init__(self, doubleClickTimeout = 200): # set timeout to zero to disable double click detection within interval
        super(GraphicsScene, self).__init__()

        self.axes = AxesGraphicsItem(None)
        self.addItem(self.axes)
        
        self.moveEventCount = 0
        self.moveItem = None
        self.moveItems = [ ]
        self.moveClickStartPos = None
        self.selectMode = False
        self.startedDrawing = False

        self.doubleClickTimer = QtCore.QTimer(self)
        self.doubleClickTimer.setInterval(doubleClickTimeout)
        self.doubleClickTimer.setSingleShot(True)
        self.doubleClickCounter = 0
        self.doubleClickTimer.timeout.connect(self.onDoubleClickTimerTimeout)
        self.doubleClickStartItem = None
        self.doubleClickStartScreenPos = (0, 0)
        self.doubleClickStartScenePos = (0, 0)

        self.lastClickModifier = None
        self.isRightClick = False

        self.dialogWidget = None

        self.axisDirection = -1
        self.pointScale = 1.0
        self.pointPixelSize = 32
        self.pointSizeFixed = True
        self.pointSceneSize = 1.0

        self.keyInfo = None

        self.matchPointsVisible = True

    def setDialogWidget(self, w):
        self.dialogWidget = w

    def getDialogWidget(self):
        return self.dialogWidget

    def warning(self, title, text = None):
        if text is None:
            text = title
        QtWidgets.QMessageBox.warning(self.dialogWidget, title, text)

    def getAxesItem(self):
        return self.axes

    def areAxesVisible(self):
        return self.axes.isVisible()

    def setAxesVisible(self, visible = True):
        self.axes.setVisible(visible)

    def getAxisDirection(self):
        return self.axisDirection

    def setAxisLeft(self, v = True):
        self.axisDirection = -1 if v else 1
        for v in self.views():
            v.updateAxisOrientation()

    def isAxisLeft(self):
        return True if self.axisDirection < 0 else False

    def _getPointTransform(self):
        s = self.pointScale if self.pointSizeFixed else self.pointSceneSize
        return QtGui.QTransform(s*self.axisDirection, 0, 0, s, 0, 0)

    def _setPointScale(self, s):
        self.pointScale = s
        self._updatePointScale()

    def _updatePointScale(self):
        t = self._getPointTransform()

        for item in self.items():
            if issubclass(type(item), PointGraphicsItemBase):
                item.setTransform(t)

    def setPointSizeFixed(self, v):
        self.pointSizeFixed = True if v else False
        self._updatePointScale()

    def isPointSizeFixed(self):
        return self.pointSizeFixed

    def setPointSizePixelSize(self, p):
        self.pointPixelSize = p

    def getPointSizePixelSize(self):
        return self.pointPixelSize

    def setPointSizeSceneSize(self, p):
        self.pointSceneSize = p

    def getPointSizeSceneSize(self):
        return self.pointSceneSize

    def onScaleChanged(self, f):
        if self.pointSizeFixed:
            self._setPointScale(self.pointPixelSize/(2.2*f))

    def mousePressEvent(self, evt):

        self.setFocus(True)

        self.lastClickModifier = self._getModifierDict(evt.modifiers())

        if evt.button() == 1:
            if evt.modifiers() & QtCore.Qt.ShiftModifier:
                self.selectMode = True
                return # Don't accept, let default handle this

            evt.accept() # note: accept to disable selection

            pos = evt.scenePos()
            self.moveEventCount = 0
            self.moveClickStartPos = pos
            self.isRightClick = False

            items = self.items(pos)
            for item in items:
                if issubclass(type(item), PointGraphicsItemBase) and item.isMovable():

                    self.moveItem = item
                    self.moveItems = [ (item, item.pos()) ]

                    for i in self.selectedItems():
                        if i is not item and issubclass(type(i), PointGraphicsItemBase) and item.isMovable():
                            self.moveItems.append( (i, i.pos()) )

                    return

            self.moveItem = None
            self.moveItems = [ ]

        elif evt.button() == 2: # Right click
            self.isRightClick = True

    def _singleClickHandlerWrapper(self, item, pos):
        from pointslayer import TriangleItem

        if item is None: # may still be a triangle
            for item in self.items(QtCore.QPointF(pos[0], pos[1])):
                if type(item) == TriangleItem:
                    self.mouseHandler_click(item, pos, self.lastClickModifier)
                    break
            else:
                self.mouseHandler_click(None, pos, self.lastClickModifier)
        else:
            self.mouseHandler_click(item, pos, self.lastClickModifier)

    def onDoubleClickTimerTimeout(self):
        if self.doubleClickCounter == 1:
            self._singleClickHandlerWrapper(self.doubleClickStartItem, self.doubleClickStartScenePos)
        elif self.doubleClickCounter == 2:
            self.mouseHandler_doubleClick(self.doubleClickStartItem, self.doubleClickStartScenePos, self.lastClickModifier)
        
        # Don't do anything if there are more clicks, just reset
        self.doubleClickCounter = 0
        self.doubleClickStartItem = None
        self.doubleClickStartScreenPos = (0, 0)
        self.doubleClickStartScenePos = (0, 0)

    def basicSingleClickHandler(self, movableItem, screenPos, scenePos):
        if not self.doubleClickTimer.isActive():
            self.doubleClickCounter = 1
            self.doubleClickStartItem = movableItem
            self.doubleClickStartScreenPos = (screenPos.x(), screenPos.y())
            self.doubleClickStartScenePos = (scenePos.x(), scenePos.y())
            self.doubleClickTimer.start()
            return

        dist2 = (screenPos.x() - self.doubleClickStartScreenPos[0])**2 + (screenPos.y() - self.doubleClickStartScreenPos[1])**2
        if dist2 <= 9: # TODO: make this more configurable
            self.doubleClickCounter += 1
            return

        # Not the same position, trigger single click, stop the timer and restart with new position
        self.doubleClickTimer.stop()
        self._singleClickHandlerWrapper(self.doubleClickStartItem, self.doubleClickStartScenePos)
        self.basicSingleClickHandler(movableItem, screenPos, scenePos)

    def mouseReleaseEvent(self, evt):

        if self.selectMode: # Don't accept, let default handle this
            self.selectMode = False # end selection
            return
        
        evt.accept()
        pos = evt.scenePos()

        if self.isRightClick:
            rightClickPoint = None
            for item in self.items(pos):
                if issubclass(type(item), PointGraphicsItemBase):
                    rightClickPoint = item
                    break
            self.mouseHandler_rightClick(rightClickPoint, (pos.x(), pos.y()), self.lastClickModifier)
            return

        if self.moveEventCount == 0: #  or self.moveItem is None:
            self.basicSingleClickHandler(self.moveItem, evt.screenPos(), pos)
        else:
            moveId = uuid.uuid4()
            for i, origPos in self.moveItems:
                itemPos = pos - self.moveClickStartPos + origPos
                scenePos = (pos.x(), pos.y()) if i == self.moveItem else None
                self.mouseHandler_moved(i, scenePos, (itemPos.x(), itemPos.y()), moveId)

            if not self.moveItems:
                self.mouseHandler_stopDrawing((pos.x(), pos.y()))

        self.moveItem = None
        self.moveItems = [ ]
        self.moveClickStartPos = None

    def mouseMoveEvent(self, evt):

        pos = evt.scenePos()
        self.onMousePosition((pos.x(), pos.y()))

        if self.selectMode:
            return

        evt.accept()

        self.moveEventCount += 1

        itemPos = None
        for i, origPos in self.moveItems:
            itemPos = pos - self.moveClickStartPos + origPos
            i.setPos(itemPos)
            i.setSelected(True)
            itemPos = (itemPos.x(), itemPos.y())

        if self.moveClickStartPos and not self.moveItems: # nothing to move, drawing mode
            if self.moveEventCount == 1: # send start drawing as well
                self.startedDrawing = self.mouseHandler_startDrawing((self.moveClickStartPos.x(), self.moveClickStartPos.y()), self.lastClickModifier)

            if self.startedDrawing:
                self.mouseHandler_drawing((pos.x(), pos.y()))

    def mouseDoubleClickEvent(self, evt):
        #print("Double click")
        evt.accept() # ignore it

    def keyPressEvent(self, evt):
        if self.focusItem():
            super(GraphicsScene, self).keyPressEvent(evt)
            return

        evt.accept()
        if not evt.isAutoRepeat():
            self.keyInfo = (evt.key(), QtGui.QKeySequence(evt.key()).toString(), evt.modifiers())

    @staticmethod
    def _getModifierDict(mod):
        modInfo = { }
        modInfo["shift"] = True if mod&QtCore.Qt.ShiftModifier else False
        modInfo["alt"] = True if mod&QtCore.Qt.AltModifier else False
        modInfo["control"] = True if mod&QtCore.Qt.ControlModifier else False
        modInfo["meta"] = True if mod&QtCore.Qt.MetaModifier else False
        return modInfo

    def keyReleaseEvent(self, evt):
        if self.focusItem():
            super(GraphicsScene, self).keyReleaseEvent(evt)
            return

        evt.accept()
        if evt.isAutoRepeat():
            return

        if self.keyInfo: # and self.keyInfo[0] == evt.key(): # this extra check fails with command key on OS X
            modInfo = self._getModifierDict(self.keyInfo[2])
            self.keyHandler_clicked(self.keyInfo[0], self.keyInfo[1], modInfo)
        self.keyInfo = None

    # Reports mouse cursor position
    def onMousePosition(self, pos):
        #print("Cursor at pos", pos)
        pass

    def onPointLabelChanged(self, layerUuid, pointUuid, oldLabel, newLabel):
        #print("onPointLabelChanged", layerUuid, pointUuid, oldLabel, newLabel)
        pass

    # Override these in derived class
    def keyHandler_clicked(self, key, keyStr, modifiers):
        #print(key, modifiers)
        pass

    def mouseHandler_moved(self, movableItem, cursorPos, itemPos, moveId):
        #print("Item", movableItem, "moved to", cursorPos, itemPos)
        pass

    def mouseHandler_click(self, movableItem, pos, modInfo):
        #print("Item", movableItem, "clicked at", pos)
        sel = movableItem.isSelected() if movableItem else False

        if sel:
            movableItem.setSelected(False)
        else:
            for i in self.selectedItems():
                i.setSelected(False)

            if movableItem:
                movableItem.setSelected(True)

    def mouseHandler_doubleClick(self, movableItem, pos, modInfo):
        #print("Item", movableItem, "double clicked at", pos)
        pass

    def mouseHandler_startDrawing(self, pos, modInfo):
        #print("START:", pos)
        return False # Indicate that we shouldn't start drawing

    def mouseHandler_drawing(self, pos):
        #print("DRAWING:", pos)
        pass

    def mouseHandler_stopDrawing(self, pos):
        #print("STOP:", pos)
        pass

    def mouseHandler_rightClick(self, pointItem, pos, modInfo):
        #print("Right click", pointItem, pos, modInfo)
        pass

    @staticmethod
    def _getWidthHeight(w, h, aspect):
        if w is None:
            if h is None:
                raise Exception("Both width and height are not set")
            w = round(h*aspect)
        else: # w is not None
            if h is None:
                h = round(w/aspect)
        return w, h

    # Note that this does not turn off any objects in the scene, and should
    # always yield an image that corresponds to the specified region as it is
    # shown (same x-axis orientation). No difference is made if the x,y
    # coordinates in bottomleft/topright are swapped
    def getSceneRegionImage(self, bottomLeft, topRight, widthPixels = None, heightPixels = None):

        center = [ (bottomLeft[i] + topRight[i])/2.0 for i in range(2) ]
        sizes = [ abs(topRight[i] - bottomLeft[i]) for i in range(2) ]
        if sizes[0] == 0 or sizes[1] == 0:
            raise Exception("Region size is zero")

        widthPixels, heightPixels = self._getWidthHeight(widthPixels, heightPixels, sizes[0]/sizes[1])
        pixmap = QtGui.QPixmap(widthPixels, heightPixels)
        pixmap.fill(QtCore.Qt.black)

        painter = QtGui.QPainter(pixmap)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)

        pixRegion = QtCore.QRectF(0, 0, widthPixels, heightPixels)
        sceneRegion = QtCore.QRectF(center[0]-sizes[0]/2.0, center[1]-sizes[1]/2.0, sizes[0], sizes[1])
        self.render(painter, pixRegion, sceneRegion)
        painter.end()

        img = pixmap.toImage()
        img = img.mirrored(self.isAxisLeft(), True)
        return img
        
    # This does take the ordering in x,y coords into account
    def getSceneRegionNumPyArray(self, bottomLeft, topRight, widthPixels = None, heightPixels = None, grayScale = False):
        import numpy as np

        img = self.getSceneRegionImage(bottomLeft, topRight, widthPixels, heightPixels)
        if grayScale:
            img = img.convertToFormat(QtGui.QImage.Format_Grayscale8)
        else:
            img = img.convertToFormat(QtGui.QImage.Format_RGB32)
        
        ptr = img.constBits()
        ptr.setsize(img.byteCount())
        if grayScale:
            arr = np.array(ptr).reshape(img.height(), img.width()).copy() # TODO: is copy necessary?
        else:
            arr = np.array(ptr).reshape(img.height(), img.width(), 4).copy() # TODO: is copy necessary?
            arr2 = np.empty((img.height(), img.width(), 3), dtype = arr.dtype)
            arr2[:,:,0] = arr[:,:,2]
            arr2[:,:,1] = arr[:,:,1]
            arr2[:,:,2] = arr[:,:,0]
            arr = arr2

        swapUD = True
        if bottomLeft[1] > topRight[1]: # need another swap
            swapUD = not swapUD

        swapLR = self.isAxisLeft()
        if bottomLeft[0] > topRight[0]:
            swapLR = not swapLR

        if swapUD:
            arr = np.flip(arr, 0)
        if swapLR:
            arr = np.flip(arr, 1)
        return arr

    def areMatchPointsVisible(self):
        return self.matchPointsVisible

    def setMatchPointsVisible(self, v):
        from imagelayer import ImageGraphicsItem

        self.matchPointsVisible = v
        for i in self.items():
            if issubclass(type(i), ImageGraphicsItem):
                i.updatePoints()

