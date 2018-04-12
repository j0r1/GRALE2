from PyQt5 import QtWidgets
from listwidgetbase import ListWidgetBase
import ui_fitslistwidget

class FITSListWidget(ListWidgetBase):
    def __init__(self, layer, item, parent = None):
        super(FITSListWidget, self).__init__(layer, item, parent)
        self.ui = ui_fitslistwidget.Ui_FITSListWidget()
        self.ui.setupUi(self)

        self.ui.m_showHideButton.clicked.connect(self.onShowHide)
        self.ui.m_visibilityBox.clicked.connect(self.onLayerVisibilityChanged)

        self.ui.m_nameEdit.textEdited.connect(self.onTextChanged)
        self.ui.m_centerRaEdit.signalNewValueEntered.connect(self.onNewRA)
        self.ui.m_centerDecEdit.signalNewValueEntered.connect(self.onNewDec)
        self.ui.m_minEdit.signalNewValueEntered.connect(self.onNewMin)
        self.ui.m_maxEdit.signalNewValueEntered.connect(self.onNewMax)

    def onShowHide(self, v):
        self.ui.m_optionsWidget.setVisible(not self.ui.m_optionsWidget.isVisible())
        self.updateSizeHint()

    def markLayerVisible(self, v):
        self.ui.m_visibilityBox.setChecked(v)

    def isLayerVisible(self):
        return self.ui.m_visibilityBox.isChecked()

    def fetchSettings(self):
        l = self.getLayer()
        ra, dec = l.getCenter()
        minValue, maxValue = l.getMinMax()

        self.ui.m_nameEdit.setText(l.getName())
        self.ui.m_centerRaEdit.setValue(float(ra))
        self.ui.m_centerDecEdit.setValue(float(dec))
        self.ui.m_minEdit.setValue(minValue)
        self.ui.m_maxEdit.setValue(maxValue)

    def onTextChanged(self, s):
        l = self.getLayer()
        self.signalLayerPropertyChanged.emit(l, "name", s)

    def onNewRA(self, v):
        l = self.getLayer()
        self.signalLayerPropertyChanged.emit(l, "centerra", v)

    def onNewDec(self, v):
        l = self.getLayer()
        self.signalLayerPropertyChanged.emit(l, "centerdec", v)

    def onNewMin(self, v):
        l = self.getLayer()
        self.signalLayerPropertyChanged.emit(l, "intensmin", v)

    def onNewMax(self, v):
        l = self.getLayer()
        self.signalLayerPropertyChanged.emit(l, "intensmax", v)

def main():
    import sys
    app = QtWidgets.QApplication(sys.argv)
    i = QtWidgets.QListWidgetItem()
    w = FITSListWidget(i)
    app.exec_()

if __name__ == "__main__":
    main()

