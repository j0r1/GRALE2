from PyQt5 import QtWidgets, QtCore, QtGui

class ListWidgetBase(QtWidgets.QFrame):
    
    s_inactiveStyleSheet = "ListWidgetBase { background-color:lightGray; }"
    s_activeStyleSheet = "ListWidgetBase { background-color:#00ff00; }"
    s_inactiveSelectedStyleSheet = "ListWidgetBase { background-color:darkGray; }"
    s_activeSelectedStyleSheet = "ListWidgetBase { background-color:#008000; }"

    signalLayerPropertyChanged = QtCore.pyqtSignal(object, str, object)
    signalVisibilityChanged = QtCore.pyqtSignal(object, bool)

    def __init__(self, layer, item, parent = None):
        super(ListWidgetBase, self).__init__(parent)
        self.item = item
        self.layer = layer
        self.isSelectedFlag = False
        self.isActiveFlag = False
        self.setActive(False)
        self.setSelected(False)

        # Use a timer to update the size hint slightly later, so that
        # the constructor of the derived class has been called completely
        t = QtCore.QTimer(self)
        t.timeout.connect(self.updateSizeHint)
        t.timeout.connect(self.fetchSettings)
        t.setInterval(0)
        t.setSingleShot(True)
        t.start()

    def getLayer(self):
        return self.layer

    def updateSizeHint(self):
        self.item.setSizeHint(self.sizeHint()) 
        self.show()

    def markLayerVisible(self, v):
        raise Exception("TODO: implement in derived class")

    def isLayerVisible(self):
        raise Exception("TODO: implement in derived class")

    def fetchSettings(self):
        raise Exception("TODO: implement in derived class")

    def onLayerVisibilityChanged(self, isVisible):
        self.signalVisibilityChanged.emit(self.layer, isVisible)

    def setActive(self, a):
        self.isActiveFlag = a
        self.updateBackground()

    def setSelected(self, s):
        self.isSelectedFlag = s
        self.updateBackground()

    def updateBackground(self):
        a, s = self.isActiveFlag, self.isSelectedFlag

        if not s:
            styleSheet = ListWidgetBase.s_activeStyleSheet if a else ListWidgetBase.s_inactiveStyleSheet
        else:
            styleSheet = ListWidgetBase.s_activeSelectedStyleSheet if a else ListWidgetBase.s_inactiveSelectedStyleSheet

        self.setStyleSheet(styleSheet)
        self.setAutoFillBackground(True)

    def isActive(self):
        return self.isActiveFlag

    def keyPressEvent(self, evt):
        evt.ignore()

    def keyReleaseEvent(self, evt):
        evt.ignore()

