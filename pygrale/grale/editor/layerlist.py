from PyQt5 import QtWidgets, QtGui, QtCore
from fitslistwidget import FITSListWidget
from rgblistwidget import RGBListWidget
from pointslistwidget import PointsListWidget
from pointslayer import PointsLayer
from imagelayer import FITSImageLayer, RGBImageLayer

class PopupMenu(QtWidgets.QMenu):
    def __init__(self, item, parent):
        super(PopupMenu, self).__init__(parent)

        self.parent = parent
        self.item = item
        deleteAction = QtWidgets.QAction("Delete layer", self)
        self.addAction(deleteAction)

        deleteAction.triggered.connect(self.onDelete)

    def onDelete(self, checked):
        self.parent._onDeleteItem(self.item)

class LayerList(QtWidgets.QListWidget):

    signalLayerDeleted = QtCore.pyqtSignal(object, int)
    signalRefreshOrder = QtCore.pyqtSignal()
    signalLayerPropertyChanged = QtCore.pyqtSignal(object, str, object)
    signalVisibilityChanged = QtCore.pyqtSignal(object, bool)

    def __init__(self, parent = None):
        super(LayerList, self).__init__(parent)

        self.itemDoubleClicked.connect(self.onItemDoubleClicked)
        self.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._onContextMenu)

    def onItemDoubleClicked(self, clickedItem):
        for i in range(self.count()):
            rowItem = self.item(i)
            w = self.itemWidget(rowItem)
            w.setActive(True) if rowItem == clickedItem else w.setActive(False)

    def fetchSettings(self):
        for i in range(self.count()):
            rowItem = self.item(i)
            w = self.itemWidget(rowItem)
            w.fetchSettings()

    def _onInternalPropertyChanged(self, layer, propertyName, propertyValue):
        self.signalLayerPropertyChanged.emit(layer, propertyName, propertyValue)

    def _onInternalVisibilityChanged(self, layer, isVisible):
        self.signalVisibilityChanged.emit(layer, isVisible)

    def addLayer(self, layer, position = -1, initialVisibility = True):
        layerToClass = { 
            PointsLayer: PointsListWidget,
            RGBImageLayer: RGBListWidget,
            FITSImageLayer: FITSListWidget
        }

        className = layerToClass[type(layer)]
        item = QtWidgets.QListWidgetItem()
        w = className(layer, item)
        w.markLayerVisible(initialVisibility)
        w.signalLayerPropertyChanged.connect(self._onInternalPropertyChanged)
        w.signalVisibilityChanged.connect(self._onInternalVisibilityChanged)

        if position < 0:
            self.addItem(item)
        else:
            self.insertItem(position, item)

        self.setItemWidget(item, w)
        return self.row(item)

    def removeLayer(self, uuid):
        for i in range(self.count()):
            it = self.item(i)
            w = self.itemWidget(self.item(i))
            if w.getLayer().getUuid() == uuid:
                r = self.row(it)
                self.takeItem(r)
                return r
        raise Exception("Specified layer ID not found")

    def _onDeleteItem(self, item):
        w = self.itemWidget(item)
        r = self.row(item)
        w.setActive(False) # Make sure it's no longer active in case of undo
        self.takeItem(r)
        self.signalLayerDeleted.emit(w.getLayer(), r)

    def _onContextMenu(self, pt):
        item = self.itemAt(pt)
        if item:
            p = PopupMenu(item, self)
            p.exec_(self.mapToGlobal(pt))

    def dropEvent(self, evt):
        super(LayerList, self).dropEvent(evt)
        self.signalRefreshOrder.emit()

    def setActiveLayer(self, uuid):
        for i in range(self.count()):
            it = self.item(i)
            w = self.itemWidget(self.item(i))
            if w.getLayer().getUuid() == uuid:
                w.setActive(True)
            else:
                w.setActive(False)

    def getActiveLayer(self):

        item = None
        for i in range(self.count()):
            it = self.item(i)
            w = self.itemWidget(self.item(i))
            if w.isActive():
                if item:
                    print("Warning: multiple active layers!")
                    w.setActive(False)
                else: 
                    item = it

        if not item:
            return None, False

        w = self.itemWidget(item)
        return w.getLayer(), w.isLayerVisible()

    def getFirstVisibleFITSItem(self):

        for i in range(self.count()):
            w = self.itemWidget(self.item(i))
            layer = w.getLayer()
            if type(layer) == FITSImageLayer and w.isLayerVisible():
                return layer
        return None

    def getLayersAndVisibilities(self):
        r = [ ]
        for i in range(self.count()):
            w = self.itemWidget(self.item(i))
            r.append( (w.getLayer(), w.isLayerVisible()) )
        return r

def main():
    import sys
    import grale.images as images
    app = QtWidgets.QApplication(sys.argv)
    w = LayerList()
    w.show()

    testFileName = "/home/jori/projects/graleeditor-hg/src/hst_12817_01_acs_wfc_f606w_drz.fits"
    fitsLayer = FITSImageLayer(testFileName)

    testFileName = "/tmp/040518_LG_dark-matter-galaxy_feat.jpg"
    #testFileName = "/home/jori/projects/a3827/eso1514a.jpg"
    rgbLayer = RGBImageLayer(testFileName)

    pointsLayer = PointsLayer()
    pointsLayer.importFromImagesData(images.ImagesData.load("/tmp/testimg.imgdata"), 0)

    w.addLayer(fitsLayer)
    w.addLayer(rgbLayer)
    w.addLayer(pointsLayer)

    app.exec_()

if __name__ == "__main__":
    main()
