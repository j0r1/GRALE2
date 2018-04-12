from PyQt5 import QtWidgets, QtGui, QtCore
import json

class DoubleLineEdit(QtWidgets.QLineEdit):

    signalNewValueEntered = QtCore.pyqtSignal(float)

    def __init__(self, parent):
        super(DoubleLineEdit, self).__init__(parent)

        self.value = 0
        self.validator = QtGui.QDoubleValidator(self)
        self.setValidator(self.validator)

        self.editingFinished.connect(self.onEnterPressed)
        self.setValue(self.value)
        self.lastValueString = self.text()

    def setRange(self, bottom, top):
        self.validator.setBottom(bottom)
        self.validator.setTop(top)

    def setValue(self, v):
        if v >=  self.validator.bottom() and v <= self.validator.top():
            self.value = v
            self.setText(json.dumps(v))
            self.lastValueString = self.text()

    def onEnterPressed(self):
        s = self.text()
        try: 
            x = float(s)
        except:
            x = 0

        #print(s, self.lastValueString)
        if self.lastValueString != s:
            self.lastValueString = s
            self.setValue(x)
            self.signalNewValueEntered.emit(x)
            #print("Emitting signalNewValueEntered",x)

    def getValue(self):
        return self.value

def main():
    import sys

    class D(DoubleLineEdit):
        def __init__(self, p):
            super(D, self).__init__(p)

            self.signalNewValueEntered.connect(self.testSlot)

        def testSlot(self, x):
            print("Received", x)

    app = QtWidgets.QApplication(sys.argv)
    d = D(None)
    d.show()
    app.exec_()

if __name__ == "__main__":
    main()
