from PyQt5 import QtWidgets
import ui_nullgriddialog

class NullGridDialog(QtWidgets.QDialog):
    def __init__(self, parent = None):
        super(NullGridDialog, self).__init__(parent)
        self.ui = ui_nullgriddialog.Ui_NullGridDialog()
        self.ui.setupUi(self)
