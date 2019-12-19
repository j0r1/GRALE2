from PyQt5 import QtWidgets, QtCore, QtGui
import ui_backprojectwidget
import pointslayer
import scenes
import base
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from grale.constants import ANGLE_ARCSEC

class BackProjectWidget(QtWidgets.QDialog):
    def __init__(self, imagePlane, parent):
        super(BackProjectWidget, self).__init__(parent)
        self.ui = ui_backprojectwidget.Ui_BackProjectWidget()
        self.ui.setupUi(self)

        self.imagePlane = imagePlane
        self.bpViews = []

        self.previousWidgetIdx = None
        self.ui.tabWidget.currentChanged.connect(self.onCurrentWidgetChanged)

        self.ui.m_pointArcsecBox.clicked.connect(self._onPointTypeChanged)
        self.ui.m_pointPixelsBox.clicked.connect(self._onPointTypeChanged)
        self.ui.m_pointArcsecSize.valueChanged.connect(self._onPointSizeArcsecChanged)
        self.ui.m_pointPixelSize.valueChanged.connect(self._onPointSizePixelChanged)
        self.ui.closeButton.clicked.connect(self.accept)

        settings = QtCore.QSettings()
        winGeom = settings.value("bpwindow/geometry")
        if winGeom is not None:
            self.restoreGeometry(winGeom)

        from tools import strToBool

        pf = settings.value("bpview/pointsizefixed")
        pf = False if not pf else strToBool(pf)
        self.ui.m_pointPixelsBox.setChecked(True) if pf else self.ui.m_pointArcsecBox.setChecked(True)
        
        psp = settings.value("bpview/pointsizepixels")
        if psp:
            self.ui.m_pointPixelSize.setValue(int(psp))
        psa = settings.value("bpview/pointsizearcsec")
        if psa:
            self.ui.m_pointArcsecSize.setValue(float(psa))
        
        print(pf, psp, psa)

    def done(self, r):
        settings = QtCore.QSettings()
        settings.setValue("bpwindow/geometry", self.saveGeometry())
        settings.setValue("bpview/pointsizefixed", self.ui.m_pointPixelsBox.isChecked())
        settings.setValue("bpview/pointsizepixels", self.ui.m_pointPixelSize.value())
        settings.setValue("bpview/pointsizearcsec", self.ui.m_pointArcsecSize.value())
        super(BackProjectWidget, self).done(r)

    def onCurrentWidgetChanged(self, idx):
        if self.previousWidgetIdx is None or self.previousWidgetIdx == idx:
            self.previousWidgetIdx = idx
            return

        # Get the center and zoom from the previous scene, and apply to the new one
        prev = self.bpViews[self.previousWidgetIdx]["view"]
        ctr, scale = prev.getCenter(), prev.getScale()

        cur = self.bpViews[idx]["view"]
        cur.centerOn(ctr)
        cur.setScale(scale)
        self.bpViews[idx]["scene"].onScaleChanged(scale)

        self.previousWidgetIdx = idx

    def addBackProjectRegion(self, bgLayer, origPointsLayer, border, srcArea):
        newPointsLayer = pointslayer.PointsLayer()

        origPts = origPointsLayer.getPoints()
        borderPoly = Polygon(border)
        for i in origPts:
            pt = origPts[i]
            xy = pt["xy"]
            if borderPoly.contains(Point(*xy)):
                bpXY = (self.imagePlane.traceThetaApproximately(np.array(xy)*ANGLE_ARCSEC)/ANGLE_ARCSEC).tolist()
                newPointsLayer.setPoint(i, bpXY, origPts[i]["label"])
                # TODO: create link between original point and backproj point

        newScene = scenes.PointsSingleLayerScene(newPointsLayer, bgLayer)
        layerItem, _ = newScene.getCurrentItemAndLayer()
        layerItem.updatePoints()
        
        newView = base.GraphicsView(newScene, parent=self.ui.tabWidget)
        self.ui.tabWidget.addTab(newView, bgLayer.getName())
        self.bpViews.append({ "view": newView, "scene": newScene})

        isPixels = self.ui.m_pointPixelsBox.isChecked()
        newScene.setPointSizePixelSize(self.ui.m_pointPixelSize.value()) if isPixels else newScene.setPointSizeSceneSize(self.ui.m_pointArcsecSize.value())
        newScene.setPointSizeFixed(self.ui.m_pointPixelsBox.isChecked())

    def setViewRange(self, bl, tr):
        if self.ui.tabWidget.count() == 0:
            return

        ctr = (bl+tr)/2
        diagArcsec = np.sum((tr-bl)**2)**0.5

        w = self.ui.tabWidget.widget(0)
        whPixels = min(w.width(),w.height())

        scale = whPixels/diagArcsec
        #print("whPixels =", whPixels, "diagArcsec =", diagArcsec, "scale=", scale)

        r = QtCore.QRectF(QtCore.QPointF(bl[0], bl[1]), QtCore.QPointF(tr[0], tr[1]))
        for i in range(self.ui.tabWidget.count()):
            w = self.ui.tabWidget.widget(i)
            s = w.scene()

            w.centerOn(ctr[0], ctr[1])
            w.setScale(scale)

            s.onScaleChanged(scale)

    def updateFromSettings(self, d):
        for x in self.bpViews:
            scene, view = x["scene"], x["view"]
            scene.setAxesVisible(d["showaxes"])
            scene.setAxisLeft(d["axisleft"])
            scene.onScaleChanged(view.getScale())

    def _onPointSizeChanged(self, isPixels, s):
        for x in self.bpViews:
            scene, view = x["scene"], x["view"]

            scene.setPointSizePixelSize(s) if isPixels else scene.setPointSizeSceneSize(s)
            scene.setPointSizeFixed(scene.isPointSizeFixed())
            scene.onScaleChanged(view.getScale())

    def _onPointSizeArcsecChanged(self, s):
        self._onPointSizeChanged(False, s)

    def _onPointSizePixelChanged(self, s):
        self._onPointSizeChanged(True, s)

    def _onPointTypeChanged(self):
        isPixels = self.ui.m_pointPixelsBox.isChecked()
        for x in self.bpViews:
            scene, view = x["scene"], x["view"]
            scene.setPointSizeFixed(isPixels)
            scene.setPointSizePixelSize(self.ui.m_pointPixelSize.value()) if isPixels else scene.setPointSizeSceneSize(self.ui.m_pointArcsecSize.value())
            scene.onScaleChanged(view.getScale())

