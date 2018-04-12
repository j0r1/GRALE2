from __future__ import print_function
from PyQt5 import QtWidgets, QtGui, QtCore
import ui_contourleveldialog

class ContourLevelDialog(QtWidgets.QDialog):
    def __init__(self, levels, pathItem, parent = None):
        super(ContourLevelDialog, self).__init__(parent)

        self.ui = ui_contourleveldialog.Ui_ContourLevelDialog()
        self.ui.setupUi(self)
        
        self.pathItem = pathItem

        minValue, maxValue = 0, len(levels)-1
        self.levels = levels
        self.ui.m_pContourLevelSlider.setMinimum(minValue)
        self.ui.m_pContourLevelSlider.setMaximum(maxValue)

        self.ui.m_pContourLevelSlider.sliderMoved.connect(self.slotSliderMoved)
        self.ui.m_pContourLevelSlider.valueChanged.connect(self.slotSliderMoved)

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
        return self.selectedContour
