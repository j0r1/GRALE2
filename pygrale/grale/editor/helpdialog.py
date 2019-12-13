from PyQt5 import QtWidgets
import ui_helpdialog

class HelpDialog(QtWidgets.QDialog):
    def __init__(self, parent, msg):
        super(HelpDialog, self).__init__(parent)

        self.ui = ui_helpdialog.Ui_HelpDialog()
        self.ui.setupUi(self)
        self.ui.m_helpEdit.setText(msg)

