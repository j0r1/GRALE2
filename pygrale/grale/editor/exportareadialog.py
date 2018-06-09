from PyQt5 import QtWidgets, QtCore
import ui_exportareadialog

class ExportAreaDialog(QtWidgets.QDialog):
    def __init__(self, areaRect, viewRect, view, areaItem, parent):
        super(ExportAreaDialog, self).__init__(parent)
        self.ui = ui_exportareadialog.Ui_ExportAreaDialog()
        self.ui.setupUi(self)

        self.areaRect = areaRect
        self.viewRect = viewRect

        self.accepted.connect(self._onAccepted)
        self.ui.m_widthPixelsBox.valueChanged.connect(self._onWidthPixelsChanged)
        self.ui.m_heightPixelsBox.valueChanged.connect(self._onHeightPixelsChanged)
        self.ui.m_widthArcsecEdit.signalNewValueEntered.connect(self._onAreaSizeChanged)
        self.ui.m_heightArcsecEdit.signalNewValueEntered.connect(self._onAreaSizeChanged)
        self.ui.m_specificAreaBox.clicked.connect(self._onNewArea)
        self.ui.m_currentViewBox.clicked.connect(self._onNewArea)
        self.ui.m_centerXEdit.signalNewValueEntered.connect(self._onNewArea)
        self.ui.m_centerYEdit.signalNewValueEntered.connect(self._onNewArea)

        from tools import strToBool

        settings = QtCore.QSettings()
        isArea = settings.value("exportdialog/exportarea")
        isArea = False if isArea is None else strToBool(isArea)
        self.ui.m_specificAreaBox.setChecked(True) if isArea else self.ui.m_currentViewBox.setChecked(True)
        
        ignoreAspect = settings.value("exportdialog/ignoreaspect")
        ignoreAspect = False if ignoreAspect is None else strToBool(ignoreAspect)
        self.ui.m_aspectBox.setChecked(False) if ignoreAspect else self.ui.m_aspectBox.setChecked(True)

        if areaRect.width() != 0 and areaRect.height() != 0:
            self.ui.m_centerXEdit.setValue((areaRect.left()+areaRect.right())/2.0)
            self.ui.m_centerYEdit.setValue((areaRect.top()+areaRect.bottom())/2.0)
            self.ui.m_widthArcsecEdit.setValue(abs(areaRect.left()-areaRect.right()))
            self.ui.m_heightArcsecEdit.setValue(abs(areaRect.top()-areaRect.bottom()))
        else:
            # Use last known settings
            centerX = settings.value("exportdialog/centerx")
            centerY = settings.value("exportdialog/centery")
            width = settings.value("exportdialog/widtharcsec")
            height = settings.value("exportdialog/heightarcsec")
            
            if centerX is None: 
                centerX = 0
            if centerY is None:
                centerY = 0
            if width is None:
                width = 10
            if height is None:
                height = 10
            self.ui.m_centerXEdit.setValue(float(centerX))
            self.ui.m_centerYEdit.setValue(float(centerY))
            self.ui.m_widthArcsecEdit.setValue(float(width))
            self.ui.m_heightArcsecEdit.setValue(float(height))

        # Note that this must be set after the export area was set
        heightPixels = settings.value("exportdialog/heightpixels")
        if heightPixels is None:
            heightPixels = 512
        self.ui.m_heightPixelsBox.setValue(int(heightPixels))
        self._onHeightPixelsChanged(int(heightPixels))

        widthPixels = settings.value("exportdialog/widthpixels")
        if widthPixels is None:
            widthPixels = 512
        self.ui.m_widthPixelsBox.setValue(int(widthPixels))
        self._onWidthPixelsChanged(int(widthPixels))

        self.areaItem = areaItem
        self.view = view
        self._onNewArea()

    def getExportRect(self):
        if self.ui.m_currentViewBox.isChecked():
            return self.viewRect

        cx = self.ui.m_centerXEdit.getValue()
        cy = self.ui.m_centerYEdit.getValue()
        w = self.ui.m_widthArcsecEdit.getValue()
        h = self.ui.m_heightArcsecEdit.getValue()
        x1, x2 = cx-w/2.0, cx+w/2.0
        y1, y2 = cy-h/2.0, cy+h/2.0

        x1, x2 = min(x1, x2), max(x1, x2)
        y1, y2 = min(y1, y2), max(y1, y2)
        return QtCore.QRectF(QtCore.QPointF(x1, y1),
                             QtCore.QPointF(x2, y2))

    def getWidthPixels(self):
        return self.ui.m_widthPixelsBox.value()

    def getHeightPixels(self):
        return self.ui.m_heightPixelsBox.value()

    def _onAccepted(self):
        settings = QtCore.QSettings()
        settings.setValue("exportdialog/exportarea", self.ui.m_specificAreaBox.isChecked())
        settings.setValue("exportdialog/ignoreaspect", not self.ui.m_aspectBox.isChecked())
        settings.setValue("exportdialog/centerx", self.ui.m_centerXEdit.getValue())
        settings.setValue("exportdialog/centery", self.ui.m_centerYEdit.getValue())
        settings.setValue("exportdialog/widtharcsec", self.ui.m_widthArcsecEdit.getValue())
        settings.setValue("exportdialog/heightarcsec", self.ui.m_heightArcsecEdit.getValue())
        settings.setValue("exportdialog/heightpixels", self.ui.m_heightPixelsBox.value())
        settings.setValue("exportdialog/widthpixels", self.ui.m_widthPixelsBox.value())

    def getAspect(self):
        if self.ui.m_currentViewBox.isChecked():
            return float(self.viewRect.width())/float(self.viewRect.height())

        return abs(self.ui.m_widthArcsecEdit.getValue()/self.ui.m_heightArcsecEdit.getValue())

    def _onWidthPixelsChanged(self, w):
        if not self.ui.m_aspectBox.isChecked():
            return

        aspect = self.getAspect()
        h = int(round(w/aspect))
        self.ui.m_heightPixelsBox.setValue(h)

    def _onHeightPixelsChanged(self, h):
        if not self.ui.m_aspectBox.isChecked():
            return
    
        aspect = self.getAspect()
        w = int(round(h*aspect))
        self.ui.m_widthPixelsBox.setValue(w)

    def _onAreaSizeChanged(self, s):
        if not self.ui.m_aspectBox.isChecked():
            return

        aspect = self.getAspect()

        if self.sender() == self.ui.m_widthArcsecEdit:
            # width changed, keep height
            h = self.ui.m_heightPixelsBox.value()
            w = int(round(h*aspect))
            self.ui.m_widthPixelsBox.setValue(w)
        else:
            w = self.ui.m_widthPixelsBox.value()
            h = int(round(w/aspect))
            self.ui.m_heightPixelsBox.setValue(h)

        self._onNewArea()

    def _onNewArea(self, checked = False):
        if not self.ui.m_specificAreaBox.isChecked():
            self.view.fitInView(self.viewRect, QtCore.Qt.KeepAspectRatio)
            self.areaItem.setVisible(False)
        else:
            x, y = self.ui.m_centerXEdit.getValue(), self.ui.m_centerYEdit.getValue()
            w, h = abs(self.ui.m_widthArcsecEdit.getValue()), abs(self.ui.m_heightArcsecEdit.getValue())
            r = QtCore.QRectF(x-w/2.0, y-h/2.0, w, h)

            self.areaItem.setRect(r)
            self.areaItem.setVisible(True)

            r.adjust(-w/15, -h/15, w/15, h/15)
            self.view.fitInView(r, QtCore.Qt.KeepAspectRatio)

