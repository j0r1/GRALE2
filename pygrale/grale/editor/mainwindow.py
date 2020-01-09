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
import uuid
from debug import log
import tools
import nullgriddialog
import exportareadialog
import grale.images as images # TODO?
from grale.constants import ANGLE_ARCSEC
import numpy as np
import pickle
import backgroundprocessdialog
import backprojretracedialog
import backprojectwidget
import backprojectsettingsdialog
import imgregionsettingsdialog

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

    def recordAddNormalPoints(self, points, addId):
        super(MultipleLayerActionStack, self).recordAddNormalPoints(points, addId)
        self.scene.updateLayerWidget()

    def recordAddTriangle(self, layerUuid, triangUuid, pointUuids, addId):
        super(MultipleLayerActionStack, self).recordAddTriangle(layerUuid, triangUuid, pointUuids, addId)
        self.scene.updateLayerWidget()

    def recordAddTriangles(self, triangles, addId):
        super(MultipleLayerActionStack, self).recordAddTriangles(triangles, addId)
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
                scene.addLayerItem(x["layer"].getUuid(), x["item"])
                scene.listWidget.addLayer(x["layer"], x["pos"], x["item"].isVisible())
            else:
                scene.removeItem(x["item"])
                scene.removeLayerItem(x["layer"].getUuid())
                pos = scene.listWidget.removeLayer(x["layer"].getUuid())
                assert(pos == x["pos"])
        elif x["cmd"] == "order":
            order = x["new"] if isRedo else x["old"]
            scene.listWidget.setOrder(order)
            scene.checkLayerOrderingAndVisibilities()
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

    def recordDeleteLayer(self, layer, layerItem, position, delUuid):
        # TODO: take delUuid into account
        self.appendToStack({
            "cmd": "dellayer",
            "layer": layer,
            "item": layerItem,
            "pos": position
        })

    def recordLayerOrderChanged(self, oldOrder, newOrder):
        self.appendToStack({
            "cmd": "order",
            "old": oldOrder,
            "new": newOrder
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

        self.checkSaveStateOnExit = True

        self.ui = ui_mainwindow.Ui_MainWindow()
        self.ui.setupUi(self)
        self.show()

        self.scene, self.view = None, None
        self._startNewFile() 

        self.ui.m_listWidget.signalLayerDeleted.connect(self._onLayerDeleted)
        self.ui.m_listWidget.signalSelectedLayersDeleted.connect(self._onSelectedLayersDeleted)
        self.ui.m_listWidget.signalRefreshOrder.connect(self._onCheckLayerOrderingAndVisibilities)
        self.ui.m_listWidget.signalLayerPropertyChanged.connect(self._onLayerPropertyChanged)
        self.ui.m_listWidget.signalVisibilityChanged.connect(self._onLayerVisibilityChanged)
        self.ui.m_listWidget.signalCenterLayerInView.connect(self._onCenterLayerInView)

        self.ui.m_addPointsButton.clicked.connect(self._addPointsLayer)
        self.ui.m_addFITSButton.clicked.connect(self._addFITSLayer)
        self.ui.m_addRGBButton.clicked.connect(self._addRGBLayer)

        self.ui.actionNew.triggered.connect(self._onNew)
        self.ui.actionLoad.triggered.connect(self._onLoad)
        self.ui.actionSave.triggered.connect(self._onSave)
        self.ui.actionSave_As.triggered.connect(self._onSaveAs)
        self.ui.actionPoints_layer_per_image.triggered.connect(self._onImportImgDat_layerPerImage)
        self.ui.actionAll_in_one_points_layer.triggered.connect(self._onImportImgDat_oneLayer)
        self.ui.actionExport_to_images_data.triggered.connect(self._onExportImgDat)
        self.ui.actionImport_JSON_file.triggered.connect(self._onImportJSON)
        self.ui.actionExport_to_JSON_file.triggered.connect(self._onExportJSON)
        self.ui.actionExport_area_view.triggered.connect(self._onExportArea)
        self.ui.actionExit.triggered.connect(self.close)
        self.ui.actionUndo.triggered.connect(self._onActionUndo)
        self.ui.actionRedo.triggered.connect(self._onActionRedo)
        self.ui.actionCopy.triggered.connect(self._onActionCopy)
        self.ui.actionPaste.triggered.connect(self._onActionPaste)
        self.ui.actionCut.triggered.connect(self._onActionCut)
        self.ui.actionCreate_null_grid.triggered.connect(self._onNullGrid)
        self.ui.actionBack_project_retrace.triggered.connect(self._onBackProjectRetrace)
        self.ui.actionBack_project_and_point_select.triggered.connect(self._onPointSelect_BackProjected)
        self.ui.actionPoint_select_no_backprojection.triggered.connect(self._onPointSelect_NoBackProject)

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

        settings = QtCore.QSettings()
        mainWinGeom = settings.value("mainwindow/geometry")
        mainWinState = settings.value("mainwindow/state")
        splitterGeom = settings.value("mainwindow/splittersizes")
        if mainWinGeom:
            self.restoreGeometry(mainWinGeom)
        if mainWinState:
            self.restoreState(mainWinState)
        if splitterGeom:
            splitterGeom = list(map(int, splitterGeom))
            self.ui.m_splitter.setSizes(splitterGeom)

        from tools import strToBool

        sa = settings.value("generalview/showaxes")
        sa = True if sa is None else strToBool(sa)
        self.ui.m_axisVisibleBox.setChecked(sa)

        al = settings.value("generalview/axesleft")
        al = True if al is None else strToBool(al)
        self.ui.m_axisLeftBox.setChecked(True) if al else self.ui.m_axisRightBox.setChecked(True)

        pf = settings.value("generalview/pointsizefixed")
        pf = False if pf is None else strToBool(pf)
        self.ui.m_pointPixelsBox.setChecked(True) if pf else self.ui.m_pointArcsecBox.setChecked(True)
        
        psp = settings.value("generalview/pointsizepixels")
        if psp:
            self.ui.m_pointPixelSize.setValue(int(psp))
        psa = settings.value("generalview/pointsizearcsec")
        if psa:
            self.ui.m_pointArcsecSize.setValue(float(psa))

        self.lastLoadedImagePlane = None

        self._onGuiSettingChanged()
        self.lastSavedState = self.getCurrentStateString()
        self.lastSaveFileName = None
        self.lastImgDatFileName = None

        timer = QtCore.QTimer(self)
        timer.setInterval(0)
        timer.setSingleShot(True)
        timer.timeout.connect(self._onStartup)
        timer.start()

    def _onStartup(self):
        # Starting up, make sure view gets keyboard
        self.view.setFocus()
        self.ui.m_zoomEdit.releaseKeyboard()
        self.ui.m_xEdit.releaseKeyboard()
        self.ui.m_yEdit.releaseKeyboard()

    def _onMousePositionChanged(self, x, y):
        self.ui.m_xEdit.setValue(x)
        self.ui.m_yEdit.setValue(y)

    def _onZoomChanged(self, s):
        self.ui.m_zoomEdit.setValue(s)

    def _onLayerDeleted(self, layer, position, delUuid = None):
        item = self.scene.getLayerItem(layer.getUuid())
        self.scene.removeItem(item)
        self.scene.getActionStack().recordDeleteLayer(layer, item, position, delUuid)
        self.scene.removeLayerItem(layer.getUuid())

    def _onSelectedLayersDeleted(self, layersAndPositions):
        # TODO: in one action stack call?
        delUuid = uuid.uuid4()
        for position, layer in layersAndPositions:
            self._onLayerDeleted(layer, position, delUuid)

    def _onCheckLayerOrderingAndVisibilities(self, oldOrder = None, newOrder = None):
        self.scene.checkLayerOrderingAndVisibilities()
        if oldOrder is not None:
            self.scene.getActionStack().recordLayerOrderChanged(oldOrder, newOrder)

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
        if curState != self.lastSavedState and self.checkSaveStateOnExit:
            if QtWidgets.QMessageBox.question(self, "Warning! Unsaved changes!", "Detected changes that were not saved, exit anyway?", 
                    defaultButton=QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.Yes:
                reallyExit = True
        else:
            reallyExit = True

        if reallyExit:
            self.view.setScene(None) # Needed to prevent crash
            evt.accept()

            settings = QtCore.QSettings()
            settings.setValue("mainwindow/geometry", self.saveGeometry())
            settings.setValue("mainwindow/state", self.saveState())
            settings.setValue("mainwindow/splittersizes", self.ui.m_splitter.sizes())
            settings.setValue("generalview/showaxes", self.ui.m_axisVisibleBox.isChecked())
            settings.setValue("generalview/axesleft", self.ui.m_axisLeftBox.isChecked())
            settings.setValue("generalview/pointsizefixed", self.ui.m_pointPixelsBox.isChecked())
            settings.setValue("generalview/pointsizepixels", self.ui.m_pointPixelSize.value())
            settings.setValue("generalview/pointsizearcsec", self.ui.m_pointArcsecSize.value())

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
        self.lastSaveFileName = None

        self.scene.signalNumberPressed.connect(self._onNumberPressed)

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

            isVisible = visibilities[i] if i < len(visibilities) else True
            self.ui.m_listWidget.markLastLayerVisible(isVisible)
            item = self.scene.getLayerItem(l.getUuid())
            item.setVisible(isVisible)

        if activeLayerUuid:
            self.ui.m_listWidget.setActiveLayer(activeLayerUuid)

        self.scene.getActionStack().clear() # Don't undo adding layers 
        self.lastSavedState = self.getCurrentStateString()
        self.lastSaveFileName = fileName

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
        if self.lastSaveFileName:
            self.saveFile(self.lastSaveFileName)
        else:
            self._onSaveAs(False)

    def status(self, s, timeout = 10):
        bar = self.statusBar()
        bar.showMessage(s, int(round(timeout*1000)))

    def saveFile(self, fileName):

        state = self.getCurrentState()
        view = self.getCurrentView()
        stateAndView = { "state": state, "view": view }
        stateAndViewStr = JSONDump(stateAndView)

        try:
            f = open(fileName, "wt")
            f.write(stateAndViewStr)
            f.close()
        except Exception as e:
            self.scene.warning("Error while saving", "Encountered a problem while saving to file '{}': {}".format(fileName, e))
            return

        self.lastSavedState = JSONDump(state)
        self.lastSaveFileName = fileName

        self.status("File '{}' saved".format(os.path.basename(fileName)))

    def _onSaveAs(self, checked):

        fileName, selectedFilter = QtWidgets.QFileDialog.getSaveFileName(self, "Specify save file name", filter="JSON files (*.json)")
        if not fileName:
            return

        self.saveFile(fileName)

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

    def importImagesData(self, imgDat, which, layerTitleFormat = "Points from images data"):
        if which == "all":
            l = pointslayer.PointsLayer()
            l.importFromImagesData(imgDat, "all")
            l.setName(layerTitleFormat)
            self.addLayer(l)
        else:
            if which >= 0:
                l = pointslayer.PointsLayer()
                l.importFromImagesData(imgDat, which)
                l.setName(layerTitleFormat)
                self.addLayer(l)
            else:
                numImages = imgDat.getNumberOfImages()
                for idx in range(numImages):
                    l = pointslayer.PointsLayer()
                    l.importFromImagesData(imgDat, idx)
                    l.setName(layerTitleFormat.format(idx))
                    self.addLayer(l)

    def importFromJSON(self, dictOrArray):

        layers = [ ]
        if "state" in dictOrArray:
            layers = dictOrArray["state"]["layers"]
        elif "layers" in dictOrArray:
            layers = dictOrArray["layers"]
        elif type(dictOrArray) == dict:
            layers = [ dictOrArray ]
        elif type(dictOrArray) == list:
            layers = dictOrArray
        else:
            raise Exception("Can't understand data to import")

        for s in layers:
            l = base.Layer.fromSettings(s)
            self.addLayer(l)

    def _onImportImgDat_layerPerImage(self, checked):
        try:
            fileName, selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self, "Select images data file", filter="Images Data files (*.imgdata *.imgdat)")
            if not fileName:
                return

            imgDat = images.ImagesData.load(fileName)
            
            shortName = os.path.basename(fileName)
            layerTitle = "Index {{}} in '{}'".format(shortName)
            
            self.importImagesData(imgDat, -1, layerTitle)
            self.lastImgDatFileName = fileName
        except Exception as e:
            self.scene.warning("Unable to import file", "Encountered a problem while importing file '{}': {}".format(fileName, e))
            return

    def _onImportImgDat_oneLayer(self, checked):
        try:
            fileName, selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self, "Select images data file", filter="Images Data files (*.imgdata *.imgdat)")
            if not fileName:
                return

            imgDat = images.ImagesData.load(fileName)
            
            shortName = os.path.basename(fileName)
            layerTitle = "All images from '{}'".format(shortName)
            
            self.importImagesData(imgDat, "all", layerTitle)
            self.lastImgDatFileName = fileName
        except Exception as e:
            self.scene.warning("Unable to import file", "Encountered a problem while importing file '{}': {}".format(fileName, e))
            return

    def _onExportImgDat(self, checked):
        pointsLeftInfo = [ ]
        try:
            splitLayers = self.ui.actionSplit_layer_into_images.isChecked()
            exportGroups = self.ui.actionExport_groups.isChecked()
            exportTimeDelays = self.ui.actionExport_time_delays.isChecked()

            layers = self._getVisiblePointsLayers()
            if not layers:
                raise Exception("No visible points layers could be detected")

            imgDat, _ = tools.layersToImagesData(layers, splitLayers, exportGroups, exportTimeDelays, pointsLeftInfo=pointsLeftInfo)
        except Exception as e:
            self.scene.warning("Error while exporting", "Encountered a problem while exporting visible points layers: {}".format(e))
            self._markRemainingPoints(pointsLeftInfo)
            return

        try:
            fileName, selectedFilter = QtWidgets.QFileDialog.getSaveFileName(self, "Specify export images data file name", self.lastImgDatFileName, filter="Images Data files (*.imgdata *.imgdat)")
            if not fileName:
                return

            imgDat.save(fileName)
            self.lastImgDatFileName = fileName
        except Exception as e:
            self.scene.warning("Error while saving", "Encountered a problem while saving to file '{}': {}".format(fileName, e))
            return

    def _markRemainingPoints(self, pointsLeftInfo):
        foundPoints = False
        for d in pointsLeftInfo:
            if d["pointsleft"]:
                foundPoints = True
                break

        if not foundPoints:
            return
       
        if not self.scene.question("Select detected remaining points?"):
            return

        self.scene.clearSelection()
        for d in pointsLeftInfo:
            layerUuid = d["layer"]
            layerItem = self.scene.getLayerItem(layerUuid)
            for pt in d["pointsleft"]:
                item = layerItem.getPointItem(pt)
                item.setSelected(True)

    def _getVisiblePointsLayers(self):
        layers = self.ui.m_listWidget.getLayersAndVisibilities()
        layers = [ l for l,v in layers if v and type(l) == pointslayer.PointsLayer ]
        return layers

    def _getVisibleLayers(self):
        layers = self.ui.m_listWidget.getLayersAndVisibilities()
        layers = [ l for l,v in layers if v ]
        return layers

    def _onImportJSON(self, checked):
        try:
            fileName, selectedFilter = QtWidgets.QFileDialog.getOpenFileName(self, "Select JSON file", filter="JSON files (*.json)")
            if not fileName:
                return

            self.importFromJSON(json.load(open(fileName, "rt")))
        except Exception as e:
            self.scene.warning("Unable to import file", "Encountered a problem while importing file '{}': {}".format(fileName, e))
            return

    def _onExportJSON(self, checked):
        try:
            layers = self._getVisibleLayers()
            if not layers:
                raise Exception("No visible layers could be detected")

            layerStr = JSONDump([ l.toSettings() for l in layers ])
        except Exception as e:
            self.scene.warning("Error while exporting", "Encountered a problem while exporting visible layers: {}".format(e))
            return

        try:
            fileName, selectedFilter = QtWidgets.QFileDialog.getSaveFileName(self, "Specify export file name", filter="JSON files (*.json)")
            if not fileName:
                return

            f = open(fileName, "wt")
            f.write(layerStr)
            f.close()
        except Exception as e:
            self.scene.warning("Error while exporting", "Encountered a problem while exporting visible layers to file '{}': {}".format(fileName, e))
            return

    def _onNullGrid(self, checked):

        rectItem = self.scene._getRectItem()
        scale = self.view.getScale()
        center = self.view.getCenter()
        
        try:
            dlg = nullgriddialog.NullGridDialog(self.view, rectItem, self)
            if not dlg.exec_():
                return

            bl, tr = dlg.getBottomLeftAndTopRight()
            bl, tr = [ bl[0]*ANGLE_ARCSEC, bl[1]*ANGLE_ARCSEC ], [ tr[0]*ANGLE_ARCSEC, tr[1]*ANGLE_ARCSEC ]
            numX, numY = dlg.getNumXYPoints()
            extraRadius = dlg.getExtraCutoutRadius()
            splitImages = dlg.getSplitImages()

            visiblePointLayers = [ l for l,v in self.ui.m_listWidget.getLayersAndVisibilities() if v and type(l) == pointslayer.PointsLayer ]
            pointsLeftInfo = [ ]
            cutoutImages = [ tools.layersToImagesData([l], splitImages, pointsLeftInfo=pointsLeftInfo)[0] for l in visiblePointLayers ]

            nullImg = images.createGridTriangles(bl, tr, numX, numY, cutoutImages, enlargeHoleOffset=extraRadius*ANGLE_ARCSEC)
            self.importImagesData(nullImg, 0, "Null space")
        except Exception as e:
            self.scene.warning("Error while creating null grid", "Encountered a problem while creating the null grid: {}".format(e))
            self._markRemainingPoints(pointsLeftInfo)
        finally:
            self.scene.removeItem(rectItem)
            self.view.setScale(scale)
            self.view.centerOn(center)


    @staticmethod
    def _getPointLayerRect(l):
        pointCoords = [ l.getPoint(p)["xy"] for p in l.getPoints() ]
        if not pointCoords:
            return QtCore.QRectF(0, 0, 0, 0)
        x = [ p[0] for p in pointCoords ]
        y = [ p[1] for p in pointCoords ]
        return QtCore.QRectF(QtCore.QPointF(min(x), min(y)), QtCore.QPointF(max(x), max(y)))

    def _getViewportRect(self):
        view = self.view
        v = view.viewport().rect()
        p1 = view.mapToScene(v.topLeft())
        p2 = view.mapToScene(v.bottomRight())
        return QtCore.QRectF(QtCore.QPointF(min(p1.x(), p2.x()), min(p1.y(), p2.y())),
                             QtCore.QPointF(max(p1.x(), p2.x()), max(p1.y(), p2.y())))

    def _onExportArea(self, checked):
        
        visiblePointLayers = [ l for l,v in self.ui.m_listWidget.getLayersAndVisibilities() if v and type(l) == pointslayer.PointsLayer ]
        r = QtCore.QRectF(0, 0, 0, 0)
        for l in visiblePointLayers:
            r |= self._getPointLayerRect(l)

        r2 = self._getViewportRect()

        try:
            rectItem = self.scene._getRectItem()
            dlg = exportareadialog.ExportAreaDialog(r, r2, self.view, rectItem, self)
            if not dlg.exec_():
                return

        finally:
            self.scene.removeItem(rectItem)

        widthPixels = dlg.getWidthPixels()
        heightPixels = dlg.getHeightPixels()

        r = dlg.getExportRect()
        centerX, centerY = (r.left() + r.right())/2.0, (r.top() + r.bottom())/2.0
        widthArcsec, heightArcsec = abs(r.left() - r.right()), abs(r.top() - r.bottom())

        try:
            img = self.scene.getSceneRegionImage([ centerX-widthArcsec/2.0, centerY-heightArcsec/2.0 ],
                                           [ centerX+widthArcsec/2.0, centerY+heightArcsec/2.0 ], widthPixels, heightPixels)
        except Exception as e:
            self.scene.warning("Error while exporting", "Encountered a problem while exporting to image: {}".format(e))
            return

        fileName, selectedFilter = QtWidgets.QFileDialog.getSaveFileName(self, "Specify export file name", filter="Images (*.jpg *.png)")
        if not fileName:
            return

        if not img.save(fileName):
            self.scene.warning("Error while exporting", "Unable to write to file {}".format(fileName))
        
    def _onCenterLayerInView(self, layer):
        pts = layer.getPoints()
        if len(pts) == 0:
            return

        pointCoords = np.array([ pts[k]["xy"] for k in pts ], dtype=np.double)
        x0, x1 = pointCoords[:,0].min(), pointCoords[:,0].max()
        y0, y1 = pointCoords[:,1].min(), pointCoords[:,1].max()
        cx = float((x0+x1)/2.0)
        cy = float((y0+y1)/2.0)
        self.view.centerOn(cx, cy)
        self.view.setFocus()

    #def keyPressEvent(self, evt):
    #    self.view.setFocus()

    def setLastImageDataFileName(self, n):
        self.lastImgDatFileName = n

    def _onNumberPressed(self, keyStr, modifiers):
        if modifiers["control"] == True and modifiers["alt"] == False:
            # Center on an image
            imageNr = int(keyStr)
            if imageNr == 0:
                imageNr = 10
            self.ui.m_listWidget.centerOnImageNumber(imageNr)

    def setCheckSaveStateOnExit(self, v):
        self.checkSaveStateOnExit = v

    def _hidePointsLayers(self, visibilities):
        for l, v in visibilities:
            if v and type(l) == pointslayer.PointsLayer:
                uuid = l.getUuid()
                layerItem = self.scene.getLayerItem(uuid)
                layerItem.setVisible(False)

    def _restorePointsLayers(self, visibilities):
        for l, v in visibilities:
            if v and type(l) == pointslayer.PointsLayer:
                uuid = l.getUuid()
                layerItem = self.scene.getLayerItem(uuid)
                layerItem.setVisible(True)

    @staticmethod
    def _createRGBLayerForImage(img, fn, layerName, bl, tr, yMirror):
        if fn is not None:
            img.save(fn)
            l = imagelayer.RGBImageLayer(fn, layerName)
        else:
            l = imagelayer.RGBImageLayer(img, layerName)

        if yMirror:
            mp = [  [ [0, img.height()], bl], 
                    [ [0, 0], [ bl[0], tr[1]] ],
                    [ [img.width(), img.height()], [tr[0], bl[1]] ],
                    [ [img.width(), 0], tr ] ]
        else:
            mp = [  [ [0, 0], bl], 
                    [ [0, img.height()], [ bl[0], tr[1]] ],
                    [ [img.width(), 0], [tr[0], bl[1]] ],
                    [ [img.width(), img.height()], tr ] ]

        l.matchToPoints(mp, True)
        return l

    def _onBackProjectRetrace(self, v):

        try:
            dlg = backprojretracedialog.BackprojRetraceDialog(self, self.lastLoadedImagePlane)
            if not dlg.exec_():
                return

            imgPlaneInfo = dlg.getImagePlaneInfo()
            if imgPlaneInfo is None:
                self.scene.warning("No image plane", "No image plane was selected")
                return

            self.lastLoadedImagePlane = imgPlaneInfo
            ip = imgPlaneInfo["imgplane"]

            splitLayers = dlg.getSplitLayersFlag()
            extra = ANGLE_ARCSEC*dlg.getExtraBorder()
            numPix = dlg.getInputImagePixels() # scaled according to aspect ratio
            numBPPix = dlg.getBackProjPixels() # same used in x and y direction no aspect ratio
            numRetracePix = dlg.getRelensPixels()  # using aspect ratio
            numResample = dlg.getResampleParameter()
            relensSeparately = dlg.getLensSeparateImages()
            overWriteFiles = dlg.getOverwriteFlag()
            bpFileNameTemplate = dlg.getBackprojectFileTemplate()
            bpLayerNameTemplate = dlg.getBackprojectLayerTemplate()
            relensFileNameTemplate = dlg.getRelensFileTemplate()
            relensLayerNameTemplate = dlg.getRelensLayerTemplate()
            newImageDir = dlg.getOutputDirectory()

            class Dlg(backgroundprocessdialog.BackgroundProcessDialog):
                def __init__(self, parent):
                    super(Dlg, self).__init__(parent, "Back-projecting and re-tracing", "Back-projecting and re-tracing...")
                    self.excepts = []
                    self.parent = parent
                    self.closed = False
                    self.newLayers = []

                def closeEvent(self, evt):
                    super(Dlg, self).closeEvent(evt)
                    self.closed = True
                    self.ui.infoLabel.setText("Cancelling...")

                def run(self):
                    def cb(msg):
                        if self.closed:
                            raise Exception("User cancelled")
                        self.ui.infoLabel.setText(msg)

                    try:
                        self.newLayers, _, _, _ = self.parent.backprojectRetrace(ip, splitLayers, extra, numPix, numBPPix, numRetracePix, numResample,
                                relensSeparately, overWriteFiles, bpFileNameTemplate, bpLayerNameTemplate,
                                relensFileNameTemplate, relensLayerNameTemplate, newImageDir, cb, addLayers = False)
                    except Exception as e:
                        self.excepts.append(str(e))

            dlg = Dlg(self)
            dlg.exec_()

            if dlg.excepts:
                raise Exception(dlg.excepts[0])

            # We need to add the layers from this thread, not from the background thread
            for l in dlg.newLayers:
                self.addLayer(l)

        except Exception as e:
            self.scene.warning("Exception occurred: {}".format(e))
            import traceback
            traceback.print_exc()
            return

    def setImagePlane(self, imgPlane, desc = ""):
        # Do some checks to verify that it's an image plane
        imgPlane.getRenderInfo()
        imgPlane.getCriticalLines()
        self.lastLoadedImagePlane = { "imgplane": imgPlane, "description": desc }

    def backprojectRetrace(self, imgPlane, splitLayers=True, extra=0, 
                           numImgPix = 1024, # Uses aspect ratio
                           numBPPix = 1024, # same used in x and y direction, is this ok? perhaps we'd lose information otherwise?
                           numRetracePix = 1024, # Again scaled according to aspect ratio
                           numResample = 1,
                           relensSeparately = True,
                           overWriteFiles = False,
                           bpFileNameTemplate = "img_{srcidx}_backproj.png",
                           bpLayerNameTemplate = "Source shape for image {srcidx}: {fn}",
                           relensFileNameTemplate = "img_{srcidx}_to_{tgtidx}_relensed.png",
                           relensLayerNameTemplate = "Relensed source from image {srcidx} to {tgtidx}: {fn}",
                           newImageDir = None,
                           progressCallback = None,
                           addLayers = True):

        ip = imgPlane
        numPix = numImgPix
        numGPUXY = numRetracePix

        if newImageDir is None:
            newImageDir = os.getcwd()
        if progressCallback is None:
            def dummyCb(msg):
                pass
            progressCallback = dummyCb

        layers = self._getVisiblePointsLayers()
        if not layers:
            raise Exception("No visible points layers could be detected")

        newLayers = [ ]
        visibilities = self.ui.m_listWidget.getLayersAndVisibilities()
        self._hidePointsLayers(visibilities)
        try:
            exportTimeDelays = False
            exportGroups = False
            pointsLeftInfo = []
            imgDat, usedLayers = tools.layersToImagesData(layers, splitLayers, exportGroups, exportTimeDelays, pointsLeftInfo=None, ignoreRemainingPoints=True)

            borders = [ ]
            for i in range(imgDat.getNumberOfImages()):
                border = None
                if imgDat.getNumberOfImagePoints(i) < 3: # Not enough points for triangulation or hull, use first point
                    border = [ imgDat.getImagePointPosition(i, 0) ]
                else:
                    for fn in [ imgDat.getBorder, imgDat.getConvexHull ]:
                        try:
                            border = fn(i)
                        except Exception as e:
                            print("Warning:", e)

                if not border:
                    raise Exception("Unable to get a border for image {}".format(i+1))
                borders.append(border)

            if not overWriteFiles: # Check what would be overwritten
                filesToWrite = [ ]
                if relensSeparately:
                    for idx in range(len(borders)):
                        filesToWrite.append(os.path.join(newImageDir, bpFileNameTemplate.format(srcidx=idx+1, tgtidx=0)))

                        for tgtidx in range(len(borders)):
                            filesToWrite.append(os.path.join(newImageDir, relensFileNameTemplate.format(srcidx=idx+1, tgtidx=tgtidx+1)))
                else:
                    for idx in range(len(borders)):
                        for tmpl in [ bpFileNameTemplate, relensFileNameTemplate ]:
                            filesToWrite.append(os.path.join(newImageDir, tmpl.format(srcidx=idx+1, tgtidx=0)))

                filesToOverwrite = [ fn for fn in filesToWrite if os.path.exists(fn) ]
                if filesToOverwrite:
                    raise Exception(f"{len(filesToOverwrite)} files would be overwritten, first are:\n" + "\n".join(filesToOverwrite[:10]))

                if len(set(filesToWrite)) != len(filesToWrite):
                    raise Exception("Some output files would overwrite each other")

            borders = [ np.array(images.enlargePolygon(b, extra))/ANGLE_ARCSEC for b in borders ]

            def addImg(img, srcidx, tgtidx, bl, tr, isSrc):
                yMirror = True if isSrc else False
                fnTemplate = bpFileNameTemplate if isSrc else relensFileNameTemplate
                layerNameTemplate = bpLayerNameTemplate if isSrc else relensLayerNameTemplate

                fn = None if fnTemplate is None else fnTemplate.format(srcidx=srcidx+1, tgtidx=tgtidx+1)
                newLayers.append(MainWindow._createRGBLayerForImage(
                    img, fn, layerNameTemplate.format(srcidx=srcidx+1, tgtidx=tgtidx+1, fn=fn),
                    bl, tr, isSrc))

            import openglhelper

            srcAreas = []
            for idx in range(len(borders)):
                border = borders[idx]
                img = self.scene.getSceneRegionImage_minMax(border.min(0), border.max(0), [ numPix, None], border)

                if ip:
                    progressCallback(f"Creating source shape from image {idx+1}")
                    centerX, centerY = 0.5*(border.max(0)+border.min(0))
                    widthArcsec, heightArcsec = border.max(0)-border.min(0)
                    imgSrc, srcbl, srctr = openglhelper.backProject(ip, img, [centerX, centerY], 
                                                      [widthArcsec, heightArcsec], [numBPPix, numBPPix])

                    srcbl, srctr = np.array(srcbl), np.array(srctr)
                else:
                    imgSrc = img.mirrored(False, True)
                    srcbl = border.min(0)
                    srctr = border.max(0)

                srcCtr = (srcbl+srctr)*0.5
                srcSize = srctr-srcbl

                srcAreas.append([srcbl, srctr])

                addImg(imgSrc, idx, -1, srcbl, srctr, True)

                if numRetracePix <= 0:
                    continue

                def getDims(tgtBl, tgtTr):
                    tgtW, tgtH = tgtTr-tgtBl 
                    return [ round(numGPUXY*tgtW/tgtH), numGPUXY ] if tgtW < tgtH else [ numGPUXY, round(numGPUXY*tgtH/tgtW) ]

                if relensSeparately:
                    for tgtidx in range(len(borders)):
                        tgtBorder = borders[tgtidx]
                        tgtBl, tgtTr = tgtBorder.min(0), tgtBorder.max(0)
    
                        dims = getDims(tgtBl, tgtTr)
                        tmpImg = imgSrc.mirrored(False, True)
                        
                        progressCallback(f"Re-tracing image {tgtidx+1} based on source shape from image {idx+1}")
                        tmpImg = openglhelper.trace(ip, tmpImg, srcCtr, srcSize, dims, numResample, tgtBl, tgtTr)
                        imgRelens = tmpImg.mirrored(False, True)

                        addImg(imgRelens, idx, tgtidx, tgtBl, tgtTr, False)
                else:

                    tgtidx = -1

                    ri = ip.getRenderInfo()
                    tgtBl = np.array(ri["bottomleft"])/ANGLE_ARCSEC
                    tgtTr = np.array(ri["topright"])/ANGLE_ARCSEC

                    dims = getDims(tgtBl, tgtTr)
                    tmpImg = imgSrc.mirrored(False, True)
                        
                    progressCallback(f"Re-tracing entire image plane based on source shape from image {idx+1}")
                    tmpImg = openglhelper.trace(ip, tmpImg, srcCtr, srcSize, dims, numResample, tgtBl, tgtTr)
                    imgRelens = tmpImg.mirrored(False, True)

                    addImg(imgRelens, idx, tgtidx, tgtBl, tgtTr, False)

        finally:
            self._restorePointsLayers(visibilities)

        # We need to do this at the end, when the original visibilities have been restored
        if addLayers:
            for l in newLayers:
                self.addLayer(l)

        return newLayers, usedLayers, borders, srcAreas

    def _onPointSelect_BackProjected(self):
        try:
            dlg = backprojectsettingsdialog.BackprojectSettingsDialog(self, self.lastLoadedImagePlane)
            if not dlg.exec_():
                return

            imgPlaneInfo = dlg.getImagePlaneInfo()
            if imgPlaneInfo is None:
                self.scene.warning("No image plane", "No image plane was selected")
                return

            self.lastLoadedImagePlane = imgPlaneInfo
            ip = imgPlaneInfo["imgplane"]

            splitLayers = dlg.getSplitLayersFlag()
            extra = ANGLE_ARCSEC*dlg.getExtraBorder()
            numPix = dlg.getInputImagePixels() # scaled according to aspect ratio
            numBPPix = dlg.getBackProjPixels() # same used in x and y direction no aspect ratio
            bpLayerNameTemplate = "Back-projected image {srcidx}"

            self._onPointSelect(ip, splitLayers, extra, numPix, numBPPix, bpLayerNameTemplate, True)

        except Exception as e:
            self.scene.warning("Exception occurred: {}".format(e))

    def _onPointSelect_NoBackProject(self):
        try:
            # TODO: change dialog
            dlg = imgregionsettingsdialog.ImageRegionSettingsDialog(self)
            if not dlg.exec_():
                return

            ip = None
            splitLayers = dlg.getSplitLayersFlag()
            extra = ANGLE_ARCSEC*dlg.getExtraBorder()
            numPix = dlg.getInputImagePixels() # scaled according to aspect ratio
            numBPPix = 1024 # is actually not used in this case
            layerNameTemplate = "Image {srcidx}"

            self._onPointSelect(ip, splitLayers, extra, numPix, numBPPix, layerNameTemplate, False)

        except Exception as e:
            self.scene.warning("Exception occurred: {}".format(e))
            import traceback
            traceback.print_exc()
            return

    def _onPointSelect(self, ip, splitLayers, extra, numPix, numBPPix, bpLayerNameTemplate, sameRegion):
        numRetracePix = -1
        numResample = -1
        relensSeparately = False
        overWriteFiles = True
        bpFileNameTemplate = None
        relensFileNameTemplate = None
        relensLayerNameTemplate = None # Note used
        newImageDir = None # Not used

        def cb(msg):
            pass

        newLayers, usedPointsLayers, borders, srcAreas = self.backprojectRetrace(ip, splitLayers, extra, numPix, numBPPix, numRetracePix, numResample,
                relensSeparately, overWriteFiles, bpFileNameTemplate, bpLayerNameTemplate,
                relensFileNameTemplate, relensLayerNameTemplate, newImageDir, cb, addLayers = False)

        maxImgs = 16
        if len(srcAreas) > maxImgs:
            raise Exception("Too many image regions ({}), limit is ({})".format(len(srcAreas), maxImgs))

        backprojectWindow = backprojectwidget.BackProjectWidget(ip, self, sameRegion)

        bl, tr = srcAreas[0]
        for l, usedPl, border, srcArea in zip(newLayers, usedPointsLayers, borders, srcAreas):
            bl, tr = np.minimum(bl, srcArea[0]), np.maximum(tr, srcArea[1])
            backprojectWindow.addBackProjectRegion(l, usedPl, border, srcArea)

        if sameRegion:
            backprojectWindow.setViewRange(bl, tr)

        backprojectWindow.updateFromSettings(self.getGlobalSettings())
        if not backprojectWindow.exec_():
            return
        
        deletedPoints = backprojectWindow.getDeletedPoints()
        deletedItems = []
        for layerUuid in deletedPoints:
            layerItem = self.scene.getLayerItem(layerUuid)
            for pt in deletedPoints[layerUuid]:
                deletedItems.append(layerItem.getPointItem(pt))

        self.scene.deleteItems(deletedItems, [], True, False)
                
        newAndChangedPoints = backprojectWindow.getImagePlanePoints()
        newPoints, changedPoints = [], []
        for layerUuid in newAndChangedPoints:
            for pt in newAndChangedPoints[layerUuid]:
                pt["layer"] = layerUuid
                pt["xy"] = pt["xy_imgplane"] # change the source plane pos by the imgplane pos
                newPoints.append(pt) if "new" in pt else changedPoints.append(pt)

        self.scene.addNewPoints(newPoints)
        self.scene.changePoints(changedPoints)


def main():
    checkQtAvailable()

    app = QtWidgets.QApplication(sys.argv)
    QtCore.QCoreApplication.setOrganizationDomain("grale2")
    QtCore.QCoreApplication.setOrganizationName("pygrale")
    QtCore.QCoreApplication.setApplicationName("editor")

    w = MainWindow()

    lastImgDataFilePrefix = "--imgdataname:"
    zoomPrefix = "--zoom:"
    nocheckPrefix = "--nocheck"
    imgplanePrefix = "--imgplane:"

    firstArg = True
    try:
        for a in sys.argv[1:]:
            if a.startswith(lastImgDataFilePrefix):
                n = a[len(lastImgDataFilePrefix):]
                w.setLastImageDataFileName(n)
            elif a.startswith(zoomPrefix):
                z = float(a[len(zoomPrefix):])
                w.setZoom(z)
            elif a == nocheckPrefix:
                w.setCheckSaveStateOnExit(False)
            elif a.startswith(imgplanePrefix):
                w.setImagePlane(pickle.load(open(a[len(imgplanePrefix):], "rb")))
            elif a.endswith(".json"):
                d = json.load(open(a, "rt"))
                loaded = False
                if type(d) == dict and ("layers" in d or "state" in d):
                    if firstArg:
                        w.loadFile(a)
                        loaded = True

                if not loaded:
                    w.importFromJSON(d)

            elif (a.endswith(".png") or a.endswith(".jpg")):
                l = imagelayer.RGBImageLayer(a)
                w.addLayer(l)

            elif a.endswith(".fits"):
                l = imagelayer.FITSImageLayer(a)
                w.addLayer(l)

            elif a == ".":
                l = pointslayer.PointsLayer()
                w.addLayer(l)

            elif (a.endswith(".imgdata") or a.endswith(".imgdat")):
                allInOne = True
                fileName = a
                if "," in a:
                    allInOne = False
                    imgIndex, fileName = a.split(",")
                    imgIndex = int(imgIndex)

                shortName = os.path.basename(fileName)
                imgDat = images.ImagesData.load(fileName)
                if allInOne:
                    which = "all"
                    layerTitle = "All images from '{}'".format(shortName)
                else:
                    which = imgIndex
                    if imgIndex >= 0:
                        layerTitle = "Index {} in '{}'".format(imgIndex, shortName)
                    else:
                        layerTitle = "Index {{}} in '{}'".format(shortName)

                w.importImagesData(imgDat, which, layerTitle)
                w.setLastImageDataFileName(fileName)
            else:
                raise Exception("Unknown argument '{}'".format(a))

            firstArg = False
    except Exception as e:
        print("Error processing arguments:", e)
        raise

    app.exec_()

if __name__ == "__main__":
    main()
