from PyQt5 import QtWidgets, QtCore, QtGui
import ui_backprojectwidget
import pointslayer
import scenes
import base
import numpy as np

class BackProjectWidget(QtWidgets.QWidget):
    def __init__(self, imagePlane):
        super(BackProjectWidget, self).__init__()
        self.ui = ui_backprojectwidget.Ui_BackProjectWidget()
        self.ui.setupUi(self)
        self.show()

        self.imagePlane = imagePlane
        self.bpViews = []

        self.previousWidgetIdx = None
        self.ui.tabWidget.currentChanged.connect(self.onCurrentWidgetChanged)

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

        self.previousWidgetIdx = idx

    def addBackProjectRegion(self, bgLayer, origPointsLayer, border, srcArea):
        newPointsLayer = pointslayer.PointsLayer()
        newScene = scenes.PointsSingleLayerScene(newPointsLayer, bgLayer)

        bl, tr = srcArea
        srcRegion = QtCore.QRectF(QtCore.QPointF(bl[0], bl[1]), QtCore.QPointF(tr[0], tr[1]))
        
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

            p = QtGui.QPen()
            p.setWidth(2)
            p.setColor(QtCore.Qt.blue)
            p.setCosmetic(True)

            ri = QtWidgets.QGraphicsRectItem(r)
            ri.setPen(p)
            s.addItem(ri)

