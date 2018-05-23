from PyQt5 import QtWidgets
import ui_exportareadialog

class ExportAreaDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super(ExportAreaDialog, self).__init__(parent)
        self.ui = ui_exportareadialog.Ui_ExportAreaDialog()
        self.ui.setupUi(self)
