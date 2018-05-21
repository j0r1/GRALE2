from PyQt5 import QtWidgets, QtCore, QtGui
import copy
import uuid

# Layer can have points (normal points, match points), possibly with other
# properties (provided by derived class) (e.g. label for point matching,
# time delay information. The 'transform' is the transformation needed to
# convert the point coordinates to scene coordinates, which is mainly for
# point matching: there, the coordinates are in pixel space and the 
# transformation converts them to scene space (where each unit is an arc second)

from cppqt.cppqt import Layer as LayerCPP

class Layer(LayerCPP):
    def __init__(self, name = "Untitled"):
        super(Layer, self).__init__(name)
    
    def createGraphicsItem(self):
        raise Exception("Implement this in derived class!")

    @staticmethod
    def fromSettings(settings):
        from pointslayer import PointsLayer
        from imagelayer import RGBImageLayer, FITSImageLayer

        if "rgbfile" in settings:
            return RGBImageLayer.fromSettings(settings)
        if "fitsfile" in settings:
            return FITSImageLayer.fromSettings(settings)
        return PointsLayer.fromSettings(settings)

from cppqt.cppqt import PointGraphicsItemBase, TriangleItem, LayerGraphicsItemBase as LayerGraphicsItemBaseCPP

class LayerGraphicsItemBase(LayerGraphicsItemBaseCPP):
    def __init__(self, layer, pointType, parent = None):
        super(LayerGraphicsItemBase, self).__init__(layer, pointType, parent)

    def arePointsVisible(self):
        s = self.scene()
        if not s:
            return True
        return s.areMatchPointsVisible()

    # TODO: rename/remove this?
    def getScenePosition(self, pos):
        p = self.getLayer().getImageTransform().map(pos[0], pos[1])
        return list(p)

    # TODO: rename/remove this?
    def getPixelPosition(self, pos):
        p = self.getLayer().getImageTransform().inverted()[0].map(pos[0], pos[1])
        return list(p)

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

from cppqt.cppqt import SceneBase

class GraphicsScene(SceneBase):

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

    def _updatePointScale(self):
        t = self._getPointTransform()
        self.setPointTransform(t);

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
            self.pointScale = self.pointPixelSize/(2.2*f)
            self._updatePointScale()

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
                item = PointGraphicsItemBase.getPointGraphicsItem(item)
                if item and item.isMovable():

                    self.moveItem = item
                    self.moveItems = [ (item, item.pos()) ]

                    if self.lastClickModifier["control"]:
                        for i in self.selectedItems():
                            i = PointGraphicsItemBase.getPointGraphicsItem(i)
                            if i and item.isMovable() and i != item:
                                self.moveItems.append( (i, i.pos()) )

                    return

            self.moveItem = None
            self.moveItems = [ ]

        elif evt.button() == 2: # Right click
            self.isRightClick = True

    def _singleClickHandlerWrapper(self, item, pos):
        if item is None: # may still be a triangle
            for item in self.items(QtCore.QPointF(pos[0], pos[1])):
                item = TriangleItem.getTriangleItem(item)
                if item:
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
        self.doubleClickStartScreenPos = [0, 0]
        self.doubleClickStartScenePos = [0, 0]

    def basicSingleClickHandler(self, movableItem, screenPos, scenePos):
        if not self.doubleClickTimer.isActive():
            self.doubleClickCounter = 1
            self.doubleClickStartItem = movableItem
            self.doubleClickStartScreenPos = [screenPos.x(), screenPos.y()]
            self.doubleClickStartScenePos = [scenePos.x(), scenePos.y()]
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
                item = PointGraphicsItemBase.getPointGraphicsItem(item)
                if item:
                    rightClickPoint = item
                    break
            self.mouseHandler_rightClick(rightClickPoint, [pos.x(), pos.y()], self.lastClickModifier)
            return

        if self.moveEventCount == 0: #  or self.moveItem is None:
            self.basicSingleClickHandler(self.moveItem, evt.screenPos(), pos)
        else:
            moveId = uuid.uuid4()
            for i, origPos in self.moveItems:
                i.setFlag(QtWidgets.QGraphicsItem.ItemSendsScenePositionChanges, False)
                itemPos = pos - self.moveClickStartPos + origPos
                scenePos = [pos.x(), pos.y()] if i == self.moveItem else None
                self.mouseHandler_moved(i, scenePos, [itemPos.x(), itemPos.y()], moveId)

            if not self.moveItems:
                self.mouseHandler_stopDrawing([pos.x(), pos.y()])

        self.moveItem = None
        self.moveItems = [ ]
        self.moveClickStartPos = None

    def mouseMoveEvent(self, evt):

        pos = evt.scenePos()
        self.onMousePosition([pos.x(), pos.y()])

        if self.selectMode:
            return

        evt.accept()

        self.moveEventCount += 1

        itemPos = None
        for i, origPos in self.moveItems:
            i.setFlag(QtWidgets.QGraphicsItem.ItemSendsScenePositionChanges)
            itemPos = pos - self.moveClickStartPos + origPos
            i.setPos(itemPos)
            i.setSelected(True)
            itemPos = [itemPos.x(), itemPos.y()]

        if self.moveClickStartPos and not self.moveItems: # nothing to move, drawing mode
            if self.moveEventCount == 1: # send start drawing as well
                self.startedDrawing = self.mouseHandler_startDrawing([self.moveClickStartPos.x(), self.moveClickStartPos.y()], self.lastClickModifier)

            if self.startedDrawing:
                self.mouseHandler_drawing([pos.x(), pos.y()])

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

    # Override these in derived class
    def keyHandler_clicked(self, key, keyStr, modifiers):
        #print(key, modifiers)
        pass

    def mouseHandler_moved(self, movableItem, cursorPos, itemPos, moveId):
        #print("Item", movableItem, "moved to", cursorPos, itemPos)
        pass

    def mouseHandler_click(self, movableItem, pos, modInfo):
        if movableItem:
            if modInfo["control"]:
                movableItem.setSelected(not movableItem.isSelected())
            else:
                sel = True if movableItem.isSelected() else False
                for i in self.selectedItems():
                    i.setSelected(False)

                if not sel:
                    movableItem.setSelected(True)

        else:
            for i in self.selectedItems():
                i.setSelected(False)

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
                i.setPointsVisible(v)

