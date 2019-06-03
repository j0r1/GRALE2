from PyQt5 import QtWidgets, QtGui, QtCore
from fitslistwidget import FITSListWidget
from rgblistwidget import RGBListWidget
from pointslistwidget import PointsListWidget
from pointslayer import PointsLayer
from imagelayer import FITSImageLayer, RGBImageLayer

class PopupMenu(QtWidgets.QMenu):
    def __init__(self, item, centerOption, delSelected, parent):
        super(PopupMenu, self).__init__(parent)

        self.parent = parent
        self.item = item

        deleteAction = QtWidgets.QAction("Delete layer", self)
        self.addAction(deleteAction)
        deleteAction.triggered.connect(self.onDelete)

        if delSelected:
            deleteSelectedAction = QtWidgets.QAction("Delete selected layers", self)
            self.addAction(deleteSelectedAction)
            deleteSelectedAction.triggered.connect(self.onDeleteSelected)

        if centerOption:
            centerAction = QtWidgets.QAction("Center in view", self)
            self.addAction(centerAction)
            centerAction.triggered.connect(self.onCenter)

    def onDelete(self, checked):
        self.parent._onDeleteItem(self.item)

    def onDeleteSelected(self, checked):
        self.parent._onDeleteSelectedItems()

    def onCenter(self, checked):
        self.parent._onCenterInView(self.item)

class LayerList(QtWidgets.QListWidget):

    signalLayerDeleted = QtCore.pyqtSignal(object, int)
    signalSelectedLayersDeleted = QtCore.pyqtSignal(object)
    signalRefreshOrder = QtCore.pyqtSignal(object, object)
    signalLayerPropertyChanged = QtCore.pyqtSignal(object, str, object)
    signalVisibilityChanged = QtCore.pyqtSignal(object, bool)
    signalCenterLayerInView = QtCore.pyqtSignal(object)

    def __init__(self, parent = None):
        super(LayerList, self).__init__(parent)

        self.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.itemDoubleClicked.connect(self.onItemDoubleClicked)
        self.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._onContextMenu)
        self.itemSelectionChanged.connect(self._onItemSelectionChanged)

        self.layerOrder = [ ]

    def onItemDoubleClicked(self, clickedItem):
        self._setActive(clickedItem)

    def _setActive(self, item):
        for i in range(self.count()):
            rowItem = self.item(i)
            w = self.itemWidget(rowItem)
            w.setActive(True) if rowItem == item else w.setActive(False)

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

        self._updateLayerOrder()
        return self.row(item)

    def _getLayerOrder(self):
        order = [ ]
        for i in range(self.count()):
            it = self.item(i)
            w = self.itemWidget(self.item(i))
            order.append(w.getLayer().getUuid())
        return order

    def _updateLayerOrder(self):
        self.layerOrder = self._getLayerOrder()

    def removeLayer(self, uuid):
        for i in range(self.count()):
            it = self.item(i)
            w = self.itemWidget(it)
            if w.getLayer().getUuid() == uuid:
                r = self.row(it)
                self.takeItem(r)
                self._updateLayerOrder()
                return r
        raise Exception("Specified layer ID not found")

    def _onDeleteItem(self, item):
        w = self.itemWidget(item)
        r = self.row(item)
        w.setActive(False) # Make sure it's no longer active in case of undo
        self.takeItem(r)
        self._updateLayerOrder()
        self.signalLayerDeleted.emit(w.getLayer(), r)

    def _onDeleteSelectedItems(self):
        layers = [ ]
        for i in range(self.count()):
            it = self.item(i)
            r = self.row(it)
            assert(r == i)

            if it.isSelected():
                w = self.itemWidget(it)
                w.setActive(False)
                l = w.getLayer()
                
                layers.append( (i, l) )

        layers = sorted(layers, key=lambda x: -x[0]) # sort from high to low layer number
        for r, l in layers:
            self.takeItem(r)

        self._updateLayerOrder()
        self.signalSelectedLayersDeleted.emit(layers)

    def _onContextMenu(self, pt):
        item = self.itemAt(pt)
        if item:
            w = self.itemWidget(item)
            isPointLayer = True if type(w) == PointsListWidget else False
            otherSelection = False
            for i in range(self.count()):
                it = self.item(i)
                w2 = self.itemWidget(it)
                if it.isSelected() and w2 is not w:
                    otherSelection = True
                    break

            p = PopupMenu(item, isPointLayer, otherSelection, self)
            p.exec_(self.mapToGlobal(pt))

    def dropEvent(self, evt):
        super(LayerList, self).dropEvent(evt)

        newOrder = self._getLayerOrder()
        assert(len(newOrder) == len(self.layerOrder))
        if newOrder != self.layerOrder:
            oldOrder = self.layerOrder
            self.layerOrder = newOrder
            self.signalRefreshOrder.emit(oldOrder, newOrder)

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

    def _onCenterInView(self, item):
        w = self.itemWidget(item)
        self._setActive(item)
        self.signalCenterLayerInView.emit(w.getLayer())

    def keyPressEvent(self, evt):
        evt.ignore()

    def keyReleaseEvent(self, evt):
        evt.ignore()

    def _onItemSelectionChanged(self):
        for i in range(self.count()):
            rowItem = self.item(i)
            w = self.itemWidget(rowItem)
            if w:
                w.setSelected(rowItem.isSelected())

    def setOrder(self, order):
        assert(len(order) == self.count())

        actLayer, actVisible = self.getActiveLayer()

        for u in range(len(order)):
            uuid = order[u]
            idx = self._find(uuid, u)
            if idx == u: # already at right position, don't do anything
                pass
            else:
                assert(idx > u)
                item = self.item(idx)
                w = self.itemWidget(item)
                layer = w.getLayer()
                vis = w.isLayerVisible()

                self.takeItem(idx)

                # Re-inserting w itself appears to cause a segfault,
                # just call addLayer again
                self.addLayer(layer, u, vis)

        self.setActiveLayer(actLayer.getUuid())

        # check again!
        for u in range(len(order)):
            w = self.itemWidget(self.item(u))
            assert(w.getLayer().getUuid() == order[u])

    def _find(self, uuid, startPos):
        for i in range(startPos, self.count()):
            w = self.itemWidget(self.item(i))
            if w.getLayer().getUuid() == uuid:
                return i
        raise Exception("Layer with uuid {} not found".format(uuid))

    def centerOnImageNumber(self, imageNumber):
        print("TODO: center on image number", imageNumber)
        idx = 0
        foundItem = None

        for i in range(self.count()):
            item = self.item(i)
            w = self.itemWidget(item)
            layer = w.getLayer()

            if type(layer) == PointsLayer:
                idx += 1
                if idx == imageNumber:
                    foundItem = item
                    break

        if foundItem is None:
            return

        self._onCenterInView(foundItem)

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
