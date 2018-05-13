from PyQt5 import QtCore, QtWidgets, QtGui
import ui_mainwindow
import sys
from checkqt import checkQtAvailable
import scenes
import base
import pointslayer
import imagelayer
import actionstack
import os
import json
from debug import log

import grale.images as images # TODO?

JSONDump = lambda s: json.dumps(s, sort_keys=True, indent=4, separators=(',', ': '))

class MultipleLayerActionStack(actionstack.ActionStack):
    def __init__(self, scene):
        super(MultipleLayerActionStack, self).__init__()
        self.scene = scene

    def recordLabelChange(self, layerUuid, pointUuid, oldLabel, newLabel):
        super(MultipleLayerActionStack, self).recordLabelChange(layerUuid, pointUuid, oldLabel, newLabel)
        self.scene.updateLayerWidget()

    #def recordPointMove(self, pointUuid, oldPosition, newPosition, moveId):
    #    super(MultipleLayerActionStack, self).recordPointMove(pointUuid, oldPosition, newPosition, moveId)

    def recordSetNormalPoint(self, layerUuid, pointUuid, oldPos, newPos, oldTd, newTd, oldLabel, newLabel):
        super(MultipleLayerActionStack, self).recordSetNormalPoint(layerUuid, pointUuid, oldPos, newPos, oldTd, newTd, oldLabel, newLabel)
        self.scene.updateLayerWidget()

    def recordAddNormalPoint(self, layerUuid, pointUuid, pos, timedelay, label, addId):
        super(MultipleLayerActionStack, self).recordAddNormalPoint(layerUuid, pointUuid, pos, timedelay, label, addId)
        self.scene.updateLayerWidget()

    def recordAddTriangle(self, layerUuid, triangUuid, pointUuids, addId):
        super(MultipleLayerActionStack, self).recordAddTriangle(layerUuid, triangUuid, pointUuids, addId)
        self.scene.updateLayerWidget()

    #def recordAddMatchPoint(self, layerUuid, pointUuid, pos, label, addId):
    #    super(MultipleLayerActionStack, self).recordAddMatchPoint(layerUuid, pointUuid, pos, label, addId)

    def recordDeletePoint(self, layerUuid, pointUuid, pos, timedelay, label, deleteId):
        super(MultipleLayerActionStack, self).recordDeletePoint(layerUuid, pointUuid, pos, timedelay, label, deleteId)
        self.scene.updateLayerWidget()

    def recordDeleteTriangle(self, layerUuid, triangUuid, pointUuids, deleteId):
        super(MultipleLayerActionStack, self).recordDeleteTriangle(layerUuid, triangUuid, pointUuids, deleteId)
        self.scene.updateLayerWidget()

    #def recordDeleteMatchPoint(self, layerUuid, pointUuid, pos, label, deleteId):
    #    super(MultipleLayerActionStack, self).recordDeleteMatchPoint(layerUuid, pointUuid, pos, label, deleteId)

    def recordCenterLocation(self, layerUuid, oldCenter, newCenter):
        super(MultipleLayerActionStack, self).recordCenterLocation(layerUuid, oldCenter, newCenter)
        self.scene.updateLayerWidget()

    def recordCenterAndMinMax(self, layerUuid, oldCenter, newCenter, oldMinMax, newMinMax):
        super(MultipleLayerActionStack, self).recordCenterAndMinMax(layerUuid, oldCenter, newCenter, oldMinMax, newMinMax)
        self.scene.updateLayerWidget()

    #def recordRGBTransform(self, layerUuid, oldTransform, newTransform):
    #    super(MultipleLayerActionStack, self).recordRGBTransform(layerUuid, oldTransform, newTransform)

    def customStackEntryUndoRedo(self, x, isRedo, scene):
        if x["cmd"] == "layername":
            item = scene.getLayerItem(x["layer"])
            layer = item.getLayer()
            name = x["new"] if isRedo else x["old"]
            layer.setName(name)
        elif (x["cmd"] == "addlayer" or x["cmd"] == "dellayer"):
            if isRedo:
                add = True if x["cmd"] == "addlayer" else False
            else:
                add = False if x["cmd"] == "addlayer" else True

            if add:
                scene.addItem(x["item"])
                scene.listWidget.addLayer(x["layer"], x["pos"], x["item"].isVisible())
            else:
                scene.removeItem(x["item"])
                pos = scene.listWidget.removeLayer(x["layer"].getUuid())
                assert(pos == x["pos"])
        else:
            raise Exception("Unknown command in '{}'".format(x))

    def recordLayerNameChange(self, layerUuid, oldName, newName):
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "layername" and prevOp["layer"] == layerUuid:
            prevOp["new"] = newName
        else:
            self.appendToStack({ 
                "cmd": "layername",
                "layer": layerUuid,
                "old": oldName,
                "new": newName
            })

    def recordAddLayer(self, layer, layerItem, position):
        self.appendToStack({
            "cmd": "addlayer",
            "layer": layer,
            "item": layerItem,
            "pos": position
        })

    def recordDeleteLayer(self, layer, layerItem, position):
        self.appendToStack({
            "cmd": "dellayer",
            "layer": layer,
            "item": layerItem,
            "pos": position
        })

class MultipleLayerScene(scenes.LayerScene):

    signalZoomChanged = QtCore.pyqtSignal(float)
    signalMousePositionChanged = QtCore.pyqtSignal(float, float)
    
    def __init__(self, listWidget, mainWidget):
        super(MultipleLayerScene, self).__init__()

        self.listWidget = listWidget
        self.mainWidget = mainWidget
        self.layerItems = { }

    def allocateActionStackInstance(self):
        return MultipleLayerActionStack(self)

    def addLayerItem(self, uuid, w):
        self.layerItems[uuid] = w

    def removeLayerItem(self, uuid):
        del self.layerItems[uuid]

    def getLayerItem(self, uuid):
        if uuid in self.layerItems:
            return self.layerItems[uuid]
        raise Exception("Specified layer {} not found".format(uuid))

    def getCurrentItemAndLayer(self):
        layer, isVisible = self.listWidget.getActiveLayer()
        if not (layer and isVisible):
            return None, None

        item = self.getLayerItem(layer.getUuid())
        if not item.isVisible():
            return None, None
        return item, layer

    def getPointMatchReferenceLayer(self):
        return self.listWidget.getFirstVisibleFITSItem()

    def onScaleChanged(self, s):
        super(MultipleLayerScene, self).onScaleChanged(s)
        self.signalZoomChanged.emit(s)

    def onMousePosition(self, xy):
        super(MultipleLayerScene, self).onMousePosition(xy)
        self.signalMousePositionChanged.emit(xy[0], xy[1])

    def undo(self):
        super(MultipleLayerScene, self).undo()
        self.listWidget.fetchSettings()

    def redo(self):
        super(MultipleLayerScene, self).redo()
        self.listWidget.fetchSettings()

    def updateLayerWidget(self):
        self.listWidget.fetchSettings()

    def checkLayerOrderingAndVisibilities(self):
        zValue = 0
        for l, v in self.listWidget.getLayersAndVisibilities():
            item = self.getLayerItem(l.getUuid())
            item.setVisible(v)
            item.setZValue(zValue)
            zValue += 1

    def getNewPointLabel(self):
        return self.mainWidget.getNewPointLabel()

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.ui = ui_mainwindow.Ui_MainWindow()
        self.ui.setupUi(self)
        self.show()

        self.scene, self.view = None, None
        self._startNewFile() 

        self.ui.m_listWidget.signalLayerDeleted.connect(self._onLayerDeleted)
        self.ui.m_listWidget.signalRefreshOrder.connect(self._onCheckLayerOrderingAndVisibilities)
        self.ui.m_listWidget.signalLayerPropertyChanged.connect(self._onLayerPropertyChanged)
        self.ui.m_listWidget.signalVisibilityChanged.connect(self._onLayerVisibilityChanged)

        self.ui.m_addPointsButton.clicked.connect(self._addPointsLayer)
        self.ui.m_addFITSButton.clicked.connect(self._addFITSLayer)
        self.ui.m_addRGBButton.clicked.connect(self._addRGBLayer)

        self.ui.actionNew.triggered.connect(self._onNew)
        self.ui.actionLoad.triggered.connect(self._onLoad)
        self.ui.actionSave.triggered.connect(self._onSave)
        self.ui.actionExit.triggered.connect(self.close)
        self.ui.actionUndo.triggered.connect(self._onActionUndo)
        self.ui.actionRedo.triggered.connect(self._onActionRedo)
        self.ui.actionCopy.triggered.connect(self._onActionCopy)
        self.ui.actionPaste.triggered.connect(self._onActionPaste)
        self.ui.actionCut.triggered.connect(self._onActionCut)

        self.ui.m_axisVisibleBox.clicked.connect(self._onGuiSettingChanged)
        self.ui.m_axisRightBox.clicked.connect(self._onGuiSettingChanged)
        self.ui.m_axisLeftBox.clicked.connect(self._onGuiSettingChanged)
        self.ui.m_pointArcsecBox.clicked.connect(self._onGuiSettingChanged)
        self.ui.m_pointPixelsBox.clicked.connect(self._onGuiSettingChanged)
        self.ui.m_pointArcsecSize.valueChanged.connect(self._onPointSizeArcsecChanged)
        self.ui.m_pointPixelSize.valueChanged.connect(self._onPointSizePixelChanged)
        self.ui.m_showMatchPoints.clicked.connect(self._onGuiSettingChanged)

        self.ui.m_zoomEdit.signalNewValueEntered.connect(self._onNewZoom)
        self.ui.m_xEdit.signalNewValueEntered.connect(self._onSetCenter)
        self.ui.m_yEdit.signalNewValueEntered.connect(self._onSetCenter)

        self.setGlobalSettings(self.getDefaultGlobalSettings())
        self.lastSavedState = self.getCurrentStateString()

    def _onMousePositionChanged(self, x, y):
        self.ui.m_xEdit.setValue(x)
        self.ui.m_yEdit.setValue(y)

    def _onZoomChanged(self, s):
        self.ui.m_zoomEdit.setValue(s)

    def _onLayerDeleted(self, layer, position):
        item = self.scene.getLayerItem(layer.getUuid())
        self.scene.removeItem(item)
        self.scene.getActionStack().recordDeleteLayer(layer, item, position)
        self.scene.removeLayerItem(layer.getUuid())

    def _onCheckLayerOrderingAndVisibilities(self):
        self.scene.checkLayerOrderingAndVisibilities()

    def _onLayerPropertyChanged(self, layer, propName, propValue):
        if propName == "name":
            oldName, newName = layer.getName(), propValue
            layer.setName(newName)
            self.scene.getActionStack().recordLayerNameChange(layer.getUuid(), oldName, newName)
        elif ( propName == "centerra" or propName == "centerdec" ): 
            oldRa, oldDec = layer.getCenter()
            newRa, newDec = propValue if propName == "centerra" else oldRa, propValue if propName == "centerdec" else oldDec
            layer.setCenter(newRa, newDec)
            self.scene.getActionStack().recordCenterLocation(layer.getUuid(), [ oldRa, oldDec ], [ newRa, newDec ])
            item = self.scene.getLayerItem(layer.getUuid())
            item.updateTransform()
        elif ( propName == "intensmin" or propName == "intensmax" ):
            oldRa, oldDec = layer.getCenter()
            oldMin, oldMax = layer.getMinMax()
            newMin, newMax = propValue if propName == "intensmin" else oldMin, propValue if propName == "intensmax" else oldMax
            layer.setMinMax(newMin, newMax)
            self.scene.getActionStack().recordCenterAndMinMax(layer.getUuid(), [ oldRa, oldDec ], [ oldRa, oldDec], 
                                                        [ oldMin, oldMax], [ newMin, newMax ])
            item = self.scene.getLayerItem(layer.getUuid())
            item.updateTransform()
            item.updatePixmap()
        else:
            raise Exception("Unknown property name '{}' in _onLayerPropertyChanged".format(propName))

        self.ui.m_listWidget.fetchSettings()

    def _onLayerVisibilityChanged(self, layer, isVisible):
        uuid = layer.getUuid()
        item = self.scene.getLayerItem(uuid)
        item.setVisible(isVisible)
    
    def _addPointsLayer(self, checked):
        self.addLayer(pointslayer.PointsLayer())

    def _addFITSLayer(self, checked):
        
        fileName, selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self, "Select FITS file", filter="FITS files (*.fits *.fit)")
        if not fileName:
            return

        try:
            l = imagelayer.FITSImageLayer(fileName, name="FITS file: {}".format(os.path.basename(fileName)))
        except Exception as e:
            self.scene.warning("Unable to use '{}' as FITS file: {}".format(fileName, e))
            return

        self.addLayer(l)

    def _addRGBLayer(self, checked):
        fileName, selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self, "Select FITS file", filter="Images (*.png *.jpg)")
        if not fileName:
            return

        try:
            l = imagelayer.RGBImageLayer(fileName, name="RGB file: {}".format(os.path.basename(fileName)))
        except Exception as e:
            self.scene.warning("Unable to use '{}' as RGB image: {}".format(fileName, e))
            return

        self.addLayer(l)

    def addLayer(self, l):
        w = l.createGraphicsItem()
        self.scene.addLayerItem(l.getUuid(), w)

        self.scene.addItem(w)
        w.updatePoints()

        pos = self.ui.m_listWidget.addLayer(l)
        self.ui.m_listWidget.setActiveLayer(l.getUuid())

        self._onCheckLayerOrderingAndVisibilities()
        self.scene.getActionStack().recordAddLayer(l, w, pos)

    def closeEvent(self, evt):
        reallyExit = False

        curState = self.getCurrentStateString()
        if curState != self.lastSavedState:
            if QtWidgets.QMessageBox.question(self, "Warning! Unsaved changes!", "Detected changes that were not saved, exit anyway?") == QtWidgets.QMessageBox.Yes:
                reallyExit = True
        else:
            reallyExit = True

        if reallyExit:
            self.view.setScene(None) # Needed to prevent crash
            evt.accept()
        else:
            evt.ignore()

    def getNewPointLabel(self):
        layer, isVisible = self.ui.m_listWidget.getActiveLayer()
        if not isVisible:
            return

        if type(layer) == imagelayer.FITSImageLayer:
            box = self.ui.m_nextFITSMatchPoint
        elif type(layer) == imagelayer.RGBImageLayer:
            box = self.ui.m_nextRGBMatchPoint
        else:
            return

        v = box.value()
        box.setValue(v+1)
        return str(v)

    def _onActionUndo(self, checked):
        self.scene.undo()

    def _onActionRedo(self, checked):
        self.scene.redo()

    def _onActionCopy(self, checked):
        self.scene.copy()

    def _onActionPaste(self, checked):
        self.scene.paste()

    def _onActionCut(self, checked):
        self.scene.cut()

    def setGlobalSettings(self, d):
        self.scene.setAxesVisible(d["showaxes"])
        self.ui.m_axisVisibleBox.setChecked(d["showaxes"])

        self.scene.setAxisLeft(d["axisleft"])
        self.ui.m_axisLeftBox.setChecked(d["axisleft"])
        self.ui.m_axisRightBox.setChecked(not d["axisleft"])

        self.scene.setPointSizePixelSize(d["pointsizepixelsize"])
        self.scene.setPointSizeSceneSize(d["pointsizearcsecsize"])
        self.scene.setPointSizeFixed(d["ispointsizefixed"])
        self.ui.m_pointPixelsBox.setChecked(d["ispointsizefixed"])
        self.ui.m_pointArcsecBox.setChecked(not d["ispointsizefixed"])
        self.ui.m_pointPixelSize.setValue(d["pointsizepixelsize"])
        self.ui.m_pointArcsecSize.setValue(d["pointsizearcsecsize"])

        self.scene.setMatchPointsVisible(d["showmatchpoints"])
        self.ui.m_showMatchPoints.setChecked(d["showmatchpoints"])

        self.ui.m_nextFITSMatchPoint.setValue(d["nextmatchpointfits"])
        self.ui.m_nextRGBMatchPoint.setValue(d["nextmatchpointrgb"])

        self.scene.onScaleChanged(self.view.getScale()) # May be needed to recalculate the point size transformation

    def getGlobalSettings(self):
        d = {
            "showaxes": self.scene.areAxesVisible(),
            "axisleft": self.scene.isAxisLeft(),
            "ispointsizefixed": self.scene.isPointSizeFixed(),
            "pointsizepixelsize": self.scene.getPointSizePixelSize(),
            "pointsizearcsecsize": self.scene.getPointSizeSceneSize(),
            "showmatchpoints": self.scene.areMatchPointsVisible(),
            "nextmatchpointrgb": self.ui.m_nextRGBMatchPoint.value(),
            "nextmatchpointfits": self.ui.m_nextFITSMatchPoint.value(),
        }
        return d

    def _onGuiSettingChanged(self, checked = None):
        d = {
            "showaxes": self.ui.m_axisVisibleBox.isChecked(),
            "axisleft": self.ui.m_axisLeftBox.isChecked(),
            "ispointsizefixed": self.ui.m_pointPixelsBox.isChecked(),
            "pointsizepixelsize": self.ui.m_pointPixelSize.value(),
            "pointsizearcsecsize": self.ui.m_pointArcsecSize.value(),
            "showmatchpoints": self.ui.m_showMatchPoints.isChecked(),
            "nextmatchpointrgb": self.ui.m_nextRGBMatchPoint.value(),
            "nextmatchpointfits": self.ui.m_nextFITSMatchPoint.value(),
        }
        self.setGlobalSettings(d)

    def _onPointSizePixelChanged(self, s):
        self.scene.setPointSizePixelSize(s)
        self.scene.setPointSizeFixed(self.scene.isPointSizeFixed()) # Causes update
        self.scene.onScaleChanged(self.view.getScale()) # May be needed to recalculate the point size transformation

    def _onPointSizeArcsecChanged(self, s):
        self.scene.setPointSizeSceneSize(s)
        self.scene.setPointSizeFixed(self.scene.isPointSizeFixed()) # Causes update
        self.scene.onScaleChanged(self.view.getScale()) # May be needed to recalculate the point size transformation

    def getDefaultGlobalSettings(self):
        d = {
            "showaxes": True,
            "axisleft": True,
            "ispointsizefixed": False,
            "pointsizepixelsize": 32,
            "pointsizearcsecsize": 0.1,
            "showmatchpoints": True,
            "nextmatchpointrgb": 0,
            "nextmatchpointfits": 0,
        }
        return d

    def _onNewZoom(self, v):
        self.setZoom(v)

    def _onSetCenter(self, dummy):
        self.setCenter(self.ui.m_xEdit.getValue(),  self.ui.m_yEdit.getValue())
    
    def setZoom(self, v):
        self.view.setScale(v)
        self.scene.onScaleChanged(v)

    def setCenter(self, x, y):
        self.view.centerOn(x, y)

    def getCurrentState(self):
        activeLayer, dummy = self.ui.m_listWidget.getActiveLayer()
        layers = [ ]
        visibilities = [ ]
        activeLayerIndex = -1

        for l, v in self.ui.m_listWidget.getLayersAndVisibilities():
            if l == activeLayer:
                activeLayerIndex = len(layers)

            layers.append(l.toSettings())
            visibilities.append(v)

        state = {
            "layers": layers,
            "visibilities": visibilities,
            "activelayerindex": activeLayerIndex,
            "globalsettings": self.getGlobalSettings(),
        }
        return state

    def getCurrentStateString(self):
        s = self.getCurrentState()
        stateStr = JSONDump(s)
        return stateStr

    def _onNew(self, checked):

        if self.getCurrentStateString() != self.lastSavedState:
            if QtWidgets.QMessageBox.question(self, "Warning! Unsaved changes!", "Detected changes that were not saved, start new anyway?") != QtWidgets.QMessageBox.Yes:
                return

        self._startNewFile()

    def _startNewFile(self):

        layerUids = [ l.getUuid() for l,v in self.ui.m_listWidget.getLayersAndVisibilities() ]
        for uuid in layerUids:
            self.ui.m_listWidget.removeLayer(uuid)

        if self.view:
            self.view.setScene(None) # Needed to prevent crash

        self.scene = MultipleLayerScene(self.ui.m_listWidget, self)
        self.scene.setDialogWidget(self)
        self.scene.signalMousePositionChanged.connect(self._onMousePositionChanged)
        self.scene.signalZoomChanged.connect(self._onZoomChanged)

        oldView = self.view if self.view else self.ui.m_dummyGraphicsView
        self.view = base.GraphicsView(self.scene, parent=self)

        self.ui.m_graphicsViewLayout.addWidget(self.view)
        self.ui.m_graphicsViewLayout.removeWidget(oldView)
        self.ui.m_zoomEdit.setValue(self.view.getScale())

        oldView.setParent(None)
        oldView.deleteLater()

        self.scene.getActionStack().clear() # Don't undo adding layers 
        self.lastSavedState = self.getCurrentStateString()

    def loadFile(self, fileName):
        self._startNewFile() # Clear everything
        stateAndView = json.load(open(fileName, "rt"))
        if not "state" in stateAndView:
            state = stateAndView
            view = None
        else:
            state = stateAndView["state"]
            view = stateAndView["view"]

        layers, visibilities = state["layers"], state["visibilities"]
        activeLayerIndex = state["activelayerindex"]
        globalSettings = state["globalsettings"]

        self.setGlobalSettings(globalSettings)

        activeLayerUuid = None
        for i in range(len(layers)):
            l = base.Layer.fromSettings(layers[i])
            if i == activeLayerIndex:
                activeLayerUuid = l.getUuid()

            self.addLayer(l)

        if activeLayerUuid:
            self.ui.m_listWidget.setActiveLayer(activeLayerUuid)

        self.scene.getActionStack().clear() # Don't undo adding layers 
        self.lastSavedState = self.getCurrentStateString()

        if view:
            self.setCurrentView(view)

    def _onLoad(self, checked):
        fileName, selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self, "Select load file name", filter="JSON files (*.json)")
        if not fileName:
            return

        if self.getCurrentStateString() != self.lastSavedState:
            if QtWidgets.QMessageBox.question(self, "Warning! Unsaved changes!", "Detected changes that were not saved, load file anyway?") != QtWidgets.QMessageBox.Yes:
                return

        try:
            self.loadFile(fileName)
        except Exception as e:
            self.scene.warning("Unable to load file", "Encountered a problem while loading file '{}': {}".format(fileName, e))
            self._startNewFile() # Don't start with a bad state
            return

    def _onSave(self, checked):
        state = self.getCurrentState()
        view = self.getCurrentView()
        stateAndView = { "state": state, "view": view }
        stateAndViewStr = JSONDump(stateAndView)

        fileName, selectedFilter = QtWidgets.QFileDialog.getSaveFileName(self, "Specify save file name", filter="JSON files (*.json)")
        if not fileName:
            return

        try:
            f = open(fileName, "wt")
            f.write(stateAndViewStr)
            f.close()
        except Exception as e:
            self.scene.warning("Error while saving", "Encountered a problem while saving to file '{}': {}".format(fileName, e))
            return

        self.lastSavedState = JSONDump(state)

    def getCurrentView(self):
        ctr = self.view.getCenter()
        scale = self.view.getScale()

        return { "center": [ ctr.x(), ctr.y() ], "zoom": scale }

    def setCurrentView(self, v):
        ctr = v["center"]
        s = v["zoom"]
        self.view.centerOn(ctr[0], ctr[1])
        self.view.setScale(s)
        self.ui.m_zoomEdit.setValue(self.view.getScale())
        self.scene.onScaleChanged(self.view.getScale()) # May be needed to recalculate the point size transformation

def main():
    checkQtAvailable()

    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()

    firstArg = True
    try:
        for a in sys.argv[1:]:
            if a.endswith(".json"):
                d = json.load(open(a, "rt"))
                if type(d) == dict and ("layers" in d or "state" in d):
                    if firstArg:
                        w.loadFile(a)
                    else:
                        for s in d["layers"]:
                            l = base.Layer.fromSettings(s)
                            w.addLayer(l)
                elif type(d) == list: # assume list of layers
                    for s in d:
                        l = base.Layer.fromSettings(s)
                        w.addLayer(l)
                else: # assume it's the settings of a single lager
                    l = base.Layer.fromSettings(d)
                    w.addLayer(l)

            elif (a.endswith(".png") or a.endswith(".jpg")):
                l = imagelayer.RGBImageLayer(a)
                w.addLayer(l)

            elif a.endswith(".fits"):
                l = imagelayer.FITSImageLayer(a)
                w.addLayer(l)

            elif a == ".":
                l = pointslayer.PointsLayer()
                w.addLayer(l)

            elif a.endswith(".imgdata"):
                allInOne = True
                fileName = a
                if "," in a:
                    allInOne = False
                    imgIndex, fileName = a.split(",")
                    imgIndex = int(imgIndex)

                shortName = os.path.basename(fileName)
                imgDat = images.ImagesData.load(fileName)
                if allInOne:
                    l = pointslayer.PointsLayer()
                    l.importFromImagesData(imgDat, "all")
                    l.setName("All images from '{}'".format(shortName))
                    w.addLayer(l)
                else:
                    if imgIndex >= 0:
                        l = pointslayer.PointsLayer()
                        l.importFromImagesData(imgDat, imgIndex)
                        l.setName("Index {} in '{}'".format(imgIndex, shortName))
                        w.addLayer(l)
                    else:
                        numImages = imgDat.getNumberOfImages()
                        for idx in range(numImages):
                            l = pointslayer.PointsLayer()
                            l.importFromImagesData(imgDat, idx)
                            l.setName("Index {} in '{}'".format(idx,shortName))
                            w.addLayer(l)
            else:
                raise Exception("Unknown argument '{}'".format(a))

            firstArg = False
    except Exception as e:
        print("Error processing arguments:", e)
        raise
        sys.exit(-1)

    app.exec_()

if __name__ == "__main__":
    main()
