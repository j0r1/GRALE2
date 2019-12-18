from PyQt5 import QtWidgets, QtCore, QtGui
import ui_backprojectwidget
import pointslayer
import scenes
import base
import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from grale.constants import ANGLE_ARCSEC

class BackProjectWidget(QtWidgets.QWidget):
    def __init__(self, imagePlane, mainWin):
        super(BackProjectWidget, self).__init__()
        self.ui = ui_backprojectwidget.Ui_BackProjectWidget()
        self.ui.setupUi(self)
        self.show()

        self.mainWin = mainWin
        self.imagePlane = imagePlane
        self.bpViews = []

        self.previousWidgetIdx = None
        self.ui.tabWidget.currentChanged.connect(self.onCurrentWidgetChanged)

    def closeEvent(self, evt):
        evt.accept()
        self.mainWin.onBackProjectWindowClosed(self)

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
            scene.setPointSizePixelSize(d["pointsizepixelsize"])
            scene.setPointSizeSceneSize(d["pointsizearcsecsize"])
            scene.setPointSizeFixed(d["ispointsizefixed"])
            scene.onScaleChanged(view.getScale())

    def onPointSizeChanged(self, isPixels, s):
        for x in self.bpViews:
            scene, view = x["scene"], x["view"]

            scene.setPointSizePixelSize(s) if isPixels else scene.setPointSizeSceneSize(s)
            scene.setPointSizeFixed(scene.isPointSizeFixed())
            scene.onScaleChanged(view.getScale())

