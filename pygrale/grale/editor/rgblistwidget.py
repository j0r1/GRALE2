from PyQt5 import QtWidgets
from listwidgetbase import ListWidgetBase
import ui_rgblistwidget

class RGBListWidget(ListWidgetBase):
    def __init__(self, layer, item, parent = None):
        super(RGBListWidget, self).__init__(layer, item, parent)
        self.ui = ui_rgblistwidget.Ui_RGBListWidget()
        self.ui.setupUi(self)

        self.ui.m_visibilityBox.clicked.connect(self.onLayerVisibilityChanged)
        self.ui.m_nameEdit.textEdited.connect(self.onTextChanged)

    def markLayerVisible(self, v):
        self.ui.m_visibilityBox.setChecked(v)

    def isLayerVisible(self):
        return self.ui.m_visibilityBox.isChecked()

    def fetchSettings(self):
        l = self.getLayer()
        self.ui.m_nameEdit.setText(l.getName())

    def onTextChanged(self, s):
        l = self.getLayer()
        self.signalLayerPropertyChanged.emit(l, "name", s)

def main():
    import sys
    app = QtWidgets.QApplication(sys.argv)
    i = QtWidgets.QListWidgetItem()
    w = RGBListWidget(i)
    app.exec_()

if __name__ == "__main__":
    main()

