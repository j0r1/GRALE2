from PyQt5 import QtWidgets, QtCore
import ui_exportareadialog

class ExportAreaDialog(QtWidgets.QDialog):
    def __init__(self, areaRect, viewRect, parent):
        super(ExportAreaDialog, self).__init__(parent)
        self.ui = ui_exportareadialog.Ui_ExportAreaDialog()
        self.ui.setupUi(self)

        self.areaRect = areaRect
        self.viewRect = viewRect

        self.ui.m_centerXEdit.setValue((areaRect.left()+areaRect.right())/2.0)
        self.ui.m_centerYEdit.setValue((areaRect.top()+areaRect.bottom())/2.0)
        self.ui.m_widthArcsecEdit.setValue(abs(areaRect.left()-areaRect.right()))
        self.ui.m_heightArcsecEdit.setValue(abs(areaRect.top()-areaRect.bottom()))

    def getExportRect(self):
        if self.ui.m_currentViewBox.isChecked():
            return self.viewRect

        cx = self.ui.m_centerXEdit.getValue()
        cy = self.ui.m_centerYEdit.getValue()
        w = self.ui.m_widthArcsecEdit.getValue()
        h = self.ui.m_heightArcsecEdit.getValue()
        x1, x2 = cx-w/2.0, cy+w/2.0
        y1, y2 = cy-h/2.0, cy+h/2.0

        x1, x2 = min(x1, x2), max(x1, x2)
        y1, y2 = min(y1, y2), max(y1, y2)
        return QtCore.QRectF(QtCore.QPointF(x1, y1),
                             QtCore.QPointF(x2, y2))

    def getWidthPixels(self):
        return self.ui.m_widthPixelsBox.value()

    def getHeightPixels(self):
        return self.ui.m_heightPixelsBox.value()

