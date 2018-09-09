from __future__ import print_function
from PyQt5 import QtWidgets, QtGui, QtCore
import grale.contourfinder as contourfinder
import ui_contourleveldialog
import numpy as np

class ContourLevelDialog(QtWidgets.QDialog):
    def __init__(self, clickPos, pathItem, rectItem, sceneRegionCallback, parent = None):
        super(ContourLevelDialog, self).__init__(parent)

        self.ui = ui_contourleveldialog.Ui_ContourLevelDialog()
        self.ui.setupUi(self)
        
        self.pathItem = pathItem
        self.rectItem = rectItem
        self.clickPos = clickPos
        self.sceneRegionCallback = sceneRegionCallback

        self.ui.m_pContourLevelSlider.sliderMoved.connect(self.slotSliderMoved)
        self.ui.m_pContourLevelSlider.valueChanged.connect(self.slotSliderMoved)

        self.ui.m_settingsButton.clicked.connect(self.onSettingsClicked)
        self.ui.m_updateButton.clicked.connect(self.updateSettings)
        self.resizeTimer = QtCore.QTimer()
        self.resizeTimer.timeout.connect(self.onResizeTimeout)
        self.resizeTimer.setSingleShot(True)
        self.resizeTimer.setInterval(0)

        self.ui.m_settingsBox.hide()
        self.resizeTimer.start()

        self.updateSettings()

    def onResizeTimeout(self):
        self.resize(self.width(), 0)

    def onSettingsClicked(self, checked = False):
        if self.ui.m_settingsBox.isVisible():
            self.ui.m_settingsBox.hide()
            self.resizeTimer.start()
        else:
            self.ui.m_settingsBox.show()

    def updateSettings(self, checked = False):

        self.levels = self._getLevels()
        minValue, maxValue = 0, len(self.levels)-1

        self.ui.m_pContourLevelSlider.setMinimum(minValue)
        self.ui.m_pContourLevelSlider.setMaximum(maxValue)

        startValue = int(round((maxValue-minValue)/3.0+minValue))
        if startValue < minValue: 
            startValue = minValue
        if startValue > maxValue:
            startValue = maxValue

        self.ui.m_pContourLevelSlider.setValue(startValue)
        self.slotSliderMoved(startValue)

    def slotSliderMoved(self, pos):
        l = self.levels[pos]
        path = QtGui.QPainterPath(QtCore.QPointF(l[0][0], l[0][1]))
        for i in range(1, len(l)):
            path.lineTo(QtCore.QPointF(l[i][0], l[i][1]))
        path.lineTo(QtCore.QPointF(l[0][0], l[0][1]))

        self.pathItem.setPath(path)
        self.selectedContour = pos

    def getSelectedContour(self):
        try:
            return self.levels[self.selectedContour]
        except Exception as e:
            print("WARNING: {}".format())
            return None

    def _getLevels(self):

        pixelSize = self.ui.m_pixelsBox.value()
        viewSize = self.ui.m_arcsecBox.value()
        pixelBlurSize = self.ui.m_blurBox.value()

        self.rectItem.setVisible(False)
        self.pathItem.setVisible(False)

        arr, bottomLeft, topRight = self.sceneRegionCallback(pixelSize, viewSize)
        
        self.rectItem.setVisible(True)
        self.pathItem.setVisible(True)

        # Blur the region to get smoother contours
        from scipy.ndimage.filters import gaussian_filter
        arr = gaussian_filter(arr, sigma=pixelBlurSize)
        
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
                if isInside(self.clickPos, part):
                    filteredContours.append(part)
                    break
    
        contours = filteredContours
        if len(contours) < 1:
            return [ ]

        self.rectItem.setRect(QtCore.QRectF(bottomLeft[0], bottomLeft[1], viewSize, viewSize))
        return contours
