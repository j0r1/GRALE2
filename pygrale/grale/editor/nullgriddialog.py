from PyQt5 import QtWidgets, QtCore
import ui_nullgriddialog

class NullGridDialog(QtWidgets.QDialog):
    def __init__(self, view, areaItem, parent = None):
        super(NullGridDialog, self).__init__(parent)
        self.ui = ui_nullgriddialog.Ui_NullGridDialog()
        self.ui.setupUi(self)

        self.accepted.connect(self._onAccepted)

        from tools import strToBool, valueFromSettings
        
        settings = QtCore.QSettings()

        self.ui.m_widthEdit.setValue(valueFromSettings(settings, "nulldialog/width", float, 100))
        self.ui.m_heightEdit.setValue(valueFromSettings(settings, "nulldialog/height", float, 100))
        self.ui.m_centerXEdit.setValue(valueFromSettings(settings, "nulldialog/centerx", float, 0))
        self.ui.m_centerYEdit.setValue(valueFromSettings(settings, "nulldialog/centery", float, 0))
        self.ui.m_numXBox.setValue(valueFromSettings(settings, "nulldialog/numx", int, 64))
        self.ui.m_numYBox.setValue(valueFromSettings(settings, "nulldialog/numy", int, 64))
        self.ui.m_cutoutEdit.setValue(valueFromSettings(settings, "nulldialog/extracutradius", float, 0.1))
        self.ui.m_splitBox.setChecked(valueFromSettings(settings, "nulldialog/splitimages", strToBool, True))

        self.ui.m_widthEdit.signalNewValueEntered.connect(self._onNewArea)
        self.ui.m_heightEdit.signalNewValueEntered.connect(self._onNewArea)
        self.ui.m_centerXEdit.signalNewValueEntered.connect(self._onNewArea)
        self.ui.m_centerYEdit.signalNewValueEntered.connect(self._onNewArea)

        self.areaItem = areaItem
        self.view = view
        self._onNewArea()

    def getBottomLeftAndTopRight(self):

        x, y = self.ui.m_centerXEdit.getValue(), self.ui.m_centerYEdit.getValue()
        w, h = abs(self.ui.m_widthEdit.getValue()), abs(self.ui.m_heightEdit.getValue())

        return [x-w/2.0, y-h/2.0], [x+w/2.0, y+h/2.0]

    def getExtraCutoutRadius(self):
        return abs(self.ui.m_cutoutEdit.getValue())
        
    def getNumXYPoints(self):
        return self.ui.m_numXBox.value(), self.ui.m_numYBox.value()

    def getSplitImages(self):
        return self.ui.m_splitBox.isChecked()

    def _onAccepted(self):

        settings = QtCore.QSettings()
        settings.setValue("nulldialog/width", self.ui.m_widthEdit.getValue())
        settings.setValue("nulldialog/height", self.ui.m_heightEdit.getValue())
        settings.setValue("nulldialog/centerx", self.ui.m_centerXEdit.getValue())
        settings.setValue("nulldialog/centery", self.ui.m_centerYEdit.getValue())
        settings.setValue("nulldialog/numx", self.ui.m_numXBox.value())
        settings.setValue("nulldialog/numy", self.ui.m_numYBox.value())
        settings.setValue("nulldialog/extracutradius", self.ui.m_cutoutEdit.getValue())
        settings.setValue("nulldialog/splitimages", self.ui.m_splitBox.isChecked())

    def _onNewArea(self):
        if not self.areaItem or not self.view:
            return

        x, y = self.ui.m_centerXEdit.getValue(), self.ui.m_centerYEdit.getValue()
        w, h = abs(self.ui.m_widthEdit.getValue()), abs(self.ui.m_heightEdit.getValue())
        r = QtCore.QRectF(x-w/2.0, y-h/2.0, w, h)

        self.areaItem.setRect(r)

        r.adjust(-w/15, -h/15, w/15, h/15)
        self.view.fitInView(r, QtCore.Qt.KeepAspectRatio)

