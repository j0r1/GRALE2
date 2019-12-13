from PyQt5 import QtWidgets, QtCore
import ui_backgroundprocessdialog
import time

class BackgroundProcessThread(QtCore.QThread):
    def __init__(self, dlg, procedure):
        super(BackgroundProcessThread, self).__init__(dlg)

        self.dlg = dlg
        self.procedure = procedure

    def run(self):
        if self.procedure is None: 
            self.dlg.run()
        else:
            self.procedure()

class BackgroundProcessDialog(QtWidgets.QDialog):

    _parentExec = QtWidgets.QDialog.exec_

    def __init__(self, parent, windowTitle = "Window Title", windowText = "Window text", procedure = None):
        super(BackgroundProcessDialog, self).__init__(parent)

        self.ui = ui_backgroundprocessdialog.Ui_BackgroundProcessDialog()
        self.ui.setupUi(self)

        self.ui.infoLabel.setText(windowText)
        self.setWindowTitle(windowTitle)

        self.procedure = procedure
        
    def closeEvent(self, e):
        e.ignore()

    # Make sure escape doesn't close the dialog
    def keyPressEvent(self, e):
        e.accept()

    def exec_(self):
        thread = BackgroundProcessThread(self, self.procedure)
        thread.finished.connect(self.accept, QtCore.Qt.QueuedConnection)
        thread.start()
        super(BackgroundProcessDialog, self).exec_()
        thread.finished.disconnect(self.accept)

    def run(self):
        print("Sleeping")
        time.sleep(10)
        print("Done")

