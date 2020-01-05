from PyQt5 import QtWidgets, QtCore, QtGui
import ui_backprojectwidget
import pointslayer
import scenes
import base
import numpy as np
import copy
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from grale.constants import ANGLE_ARCSEC

class BackProjectWidget(QtWidgets.QDialog):
    def __init__(self, imagePlane, parent, sameRegion):
        super(BackProjectWidget, self).__init__(parent)
        self.ui = ui_backprojectwidget.Ui_BackProjectWidget()
        self.ui.setupUi(self)

        self.imagePlane = imagePlane
        self.sameRegion = sameRegion
        self.bpViews = []
        self.originalPointInfo = { }
        self.imgplanePoints = { }
        self.deletedPoints = { }
        self.originalImagePlanePositions = { }

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

    def done(self, r):
        settings = QtCore.QSettings()
        settings.setValue("bpwindow/geometry", self.saveGeometry())
        settings.setValue("bpview/pointsizefixed", self.ui.m_pointPixelsBox.isChecked())
        settings.setValue("bpview/pointsizepixels", self.ui.m_pointPixelSize.value())
        settings.setValue("bpview/pointsizearcsec", self.ui.m_pointArcsecSize.value())
        
        if r:
            if not self._checkPoints():
                return
        super(BackProjectWidget, self).done(r)

    def onCurrentWidgetChanged(self, idx):
        if self.previousWidgetIdx is None or self.previousWidgetIdx == idx:
            self.previousWidgetIdx = idx
            return

        # Get the center and zoom from the previous scene, and apply to the new one
        if self.sameRegion:
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
            self.originalImagePlanePositions[i] = { "xy": [ xy[0], xy[1] ], "label": pt["label"], "timedelay": pt["timedelay"] }

            if borderPoly.contains(Point(*xy)):
                bpXY = (self.imagePlane.traceThetaApproximately(np.array(xy)*ANGLE_ARCSEC)/ANGLE_ARCSEC).tolist()
                newPointsLayer.setPoint(i, bpXY, pt["label"])
        
                origLayerId = origPointsLayer.getUuid()
                if not origLayerId in self.originalPointInfo:
                    self.originalPointInfo[origLayerId] = { }
                self.originalPointInfo[origLayerId][i] = newPointsLayer.getPoint(i) # Use this to check what to update

        newScene = scenes.PointsSingleLayerScene(newPointsLayer, bgLayer)
        layerItem, _ = newScene.getCurrentItemAndLayer()
        layerItem.updatePoints()
        
        newView = base.GraphicsView(newScene, parent=self.ui.tabWidget)
        self.ui.tabWidget.addTab(newView, bgLayer.getName())
        self.bpViews.append({ "view": newView, "scene": newScene, "newlayer": newPointsLayer, 
                              "origlayer": origPointsLayer, "srcarea": srcArea, "border": border })

        isPixels = self.ui.m_pointPixelsBox.isChecked()
        newScene.setPointSizePixelSize(self.ui.m_pointPixelSize.value()) if isPixels else newScene.setPointSizeSceneSize(self.ui.m_pointArcsecSize.value())
        newScene.setPointSizeFixed(self.ui.m_pointPixelsBox.isChecked())

        self._setViewRange(newView, srcArea[0], srcArea[1])

    def _checkPoints(self):
        try:
            pointsOutsideSrcArea = [ ]
            pointsOutsideImgBorder = [ ]
            pointsWithTraceException = [ ]

            pointsPerLayer = { }
            allPointsPerLayer = { }
            deletedPoints = { }

            for x in self.bpViews:
                layer, origLayer, srcArea, border = x["newlayer"], x["origlayer"], x["srcarea"], x["border"]
                bl, tr = srcArea
                borderPoly = Polygon(border)

                origLayerUuid = origLayer.getUuid()
                if not origLayerUuid in pointsPerLayer:
                    pointsPerLayer[origLayerUuid] = [ ]
                if not origLayerUuid in allPointsPerLayer:
                    allPointsPerLayer[origLayerUuid] = set()

                points = layer.getPoints()
                for uuid in points:
                    allPointsPerLayer[origLayerUuid].add(uuid)
                    pt = copy.copy(points[uuid])

                    if uuid in self.originalPointInfo[origLayerUuid]:
                        if pt == self.originalPointInfo[origLayerUuid][uuid]: 
                            # nothing changed, ignore this
                            #print("Nothing changed for", pt)
                            continue
                    else:
                        pt["new"] = True

                    pt["point"] = uuid
                    xy = pt["xy"]

                    if not (xy[0] >= bl[0] and xy[0] <= tr[0] and xy[1] >= bl[1] and xy[1] <= tr[1]):
                        pointsOutsideSrcArea.append({"point": uuid, "newlayer": layer.getUuid()})
                    else:
                        # Check if src plane positions are still the same
                        if not "new" in pt and xy == self.originalPointInfo[origLayerUuid][uuid]["xy"]:
                            # if so, just use the stored image plane pos
                            #print("Source plane pos unchanged, using image plane pos")
                            pt["xy_imgplane"] = self.originalImagePlanePositions[uuid]["xy"]
                            pointsPerLayer[origLayerUuid].append(pt)
                        else:
                            try:
                                # Trace it
                                #print("Tracing src plane pos for", uuid, "pos", xy)
                                thetas = self.imagePlane.traceBeta(np.array(xy)*ANGLE_ARCSEC)
                                #print("Thetas", thetas)
                                for theta in thetas:
                                    theta /= ANGLE_ARCSEC
                                    if borderPoly.contains(Point(theta.tolist())): # Ok, this is the one
                                        pt["xy_imgplane"] = theta
                                        pointsPerLayer[origLayerUuid].append(pt)
                                        break

                                else:
                                    pointsOutsideImgBorder.append({"point": uuid, "newlayer": layer.getUuid()})

                            except Exception as e:
                                pointsWithTraceException.append({"point": uuid, "newlayer": layer.getUuid()})
                                print(e)

            for origLayerUuid in self.originalPointInfo:
                if not origLayerUuid in allPointsPerLayer: # all points in layer are deleted
                    deletedPoints[origLayerUuid] = self.originalPointInfo[origLayerUuid]
                else:
                    allPts = allPointsPerLayer[origLayerUuid]
                    origPts = self.originalPointInfo[origLayerUuid]
                    for uuid in origPts:
                        if not uuid in allPts: # this point was removed from the layer

                            if not origLayerUuid in deletedPoints:
                                deletedPoints[origLayerUuid] = [ ]
                            deletedPoints[origLayerUuid].append(uuid)
                    
            if False:
                import pprint
                print("Points with trace exception:")
                pprint.pprint(pointsWithTraceException)
                print("Points outside image border:")
                pprint.pprint(pointsOutsideImgBorder)
                print("Points outside src area:")
                pprint.pprint(pointsOutsideSrcArea)
                print("Points per layer:")
                pprint.pprint(pointsPerLayer)
                print("Deleted points:")
                pprint.pprint(deletedPoints)

            # TODO: better messages? mark points that are problematic
            if pointsWithTraceException:
                raise Exception("A number of points ({}) could not be traced back to the image plane".format(len(pointsWithTraceException)))
            if pointsOutsideSrcArea:
                raise Exception("A number of points ({}) lie outside the back-projected region".format(len(pointsOutsideSrcArea)))
            if pointsOutsideImgBorder:
                raise Exception("A number of re-traced points ({}) couldn't be found inside the image region".format(len(pointsOutsideImgBorder)))

            for lUuid in pointsPerLayer:
                for pt in pointsPerLayer[lUuid]:
                    if "new" in pt:
                        continue

                    orig = self.originalImagePlanePositions[pt["point"]]
                    pt["xy_orig"] = orig["xy"]
                    pt["label_orig"] = orig["label"]
                    pt["timedelay_orig"] = orig["timedelay"]

                    for k in [ "xy", "xy_imgplane", "xy_orig" ]:
                        if type(pt[k]) != list:
                            pt[k] = pt[k].tolist()
        
            self.deletedPoints = deletedPoints
            self.imgplanePoints = pointsPerLayer

        except Exception as e:
            QtWidgets.QMessageBox.warning(self, "Error checking points", str(e))
            return False

        return True

    def _setViewRange(self, w, bl, tr):
        ctr = (bl+tr)/2
        diagArcsec = np.sum((tr-bl)**2)**0.5

        whPixels = min(w.width(),w.height())
        scale = whPixels/diagArcsec

        s = w.scene()
        w.centerOn(ctr[0], ctr[1])
        w.setScale(scale)

        s.onScaleChanged(scale)

    def setViewRange(self, bl, tr):
        if self.ui.tabWidget.count() == 0:
            return

        for i in range(self.ui.tabWidget.count()):
            w = self.ui.tabWidget.widget(i)
            self._setViewRange(w, bl, tr)

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
        psa = self.ui.m_pointArcsecSize.value()
        psp = self.ui.m_pointPixelSize.value()
        print(isPixels, psp, psa)
        for x in self.bpViews:
            scene, view = x["scene"], x["view"]
            scene.setPointSizePixelSize(psp) if isPixels else scene.setPointSizeSceneSize(psa)
            scene.setPointSizeFixed(isPixels)
            scene.onScaleChanged(view.getScale())

    def getImagePlanePoints(self):
        return self.imgplanePoints

    def getDeletedPoints(self):
        return self.deletedPoints

