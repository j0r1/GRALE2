from PyQt5 import QtWidgets
from listwidgetbase import ListWidgetBase
import ui_pointslistwidget

class PointsListWidget(ListWidgetBase):
    def __init__(self, layer, item, parent = None):
        super(PointsListWidget, self).__init__(layer, item, parent)
        self.ui = ui_pointslistwidget.Ui_PointsListWidget()
        self.ui.setupUi(self)

        self.ui.m_showHideButton.clicked.connect(self.onShowHide)
        self.ui.m_visibilityBox.clicked.connect(self.onLayerVisibilityChanged)
        self.ui.m_nameEdit.textEdited.connect(self.onTextChanged)

    def onShowHide(self, v):
        self.ui.m_optionsWidget.setVisible(not self.ui.m_optionsWidget.isVisible())
        self.updateSizeHint()

    def markLayerVisible(self, v):
        self.ui.m_visibilityBox.setChecked(v)

    def isLayerVisible(self):
        return self.ui.m_visibilityBox.isChecked()

    def fetchSettings(self):
        l = self.getLayer()
        numPoints = len(l.getPoints())
        numTriangles = len(l.getTriangles())

        self.ui.m_nameEdit.setText(l.getName())
        self.ui.m_numPointsEdit.setText(str(numPoints))
        self.ui.m_numTrianglesEdit.setText(str(numTriangles))

    def onTextChanged(self, s):
        l = self.getLayer()
        self.signalLayerPropertyChanged.emit(l, "name", s)

def main():
    import sys
    app = QtWidgets.QApplication(sys.argv)
    i = QtWidgets.QListWidgetItem()
    w = PointsListWidget(i)
    app.exec_()

if __name__ == "__main__":
    main()

