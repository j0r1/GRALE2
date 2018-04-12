from PyQt5 import QtWidgets
import ui_pointinfodialog

class PointInfoDialog(QtWidgets.QDialog):
    def __init__(self, pointInfo, parent = None):
        super(PointInfoDialog, self).__init__(parent)

        self.ui = ui_pointinfodialog.Ui_PointInfoDialog()
        self.ui.setupUi(self)

        self.ui.m_xEdit.setValue(pointInfo["xy"][0])
        self.ui.m_yEdit.setValue(pointInfo["xy"][1])
        if pointInfo["label"]:
            self.ui.m_groupName.setText(pointInfo["label"])
        if pointInfo["timedelay"] is not None:
            self.ui.m_tdCheckBox.setChecked(True)
            self.ui.m_tdEdit.setValue(pointInfo["timedelay"])

    def getPointInfo(self):
        ptInfo = {
            "xy": [ self.ui.m_xEdit.getValue(), self.ui.m_yEdit.getValue() ],
            "label": None if not self.ui.m_groupName.text().strip() else self.ui.m_groupName.text(),
            "timedelay": None if not self.ui.m_tdCheckBox.isChecked() else self.ui.m_tdEdit.getValue()
        }
        return ptInfo

def main():
    import sys
    import pprint

    app = QtWidgets.QApplication(sys.argv)
    ptInfo = { "xy": [1,2], "label": "hey", "timedelay": 123 }
    w = PointInfoDialog(ptInfo)
    w.show()
    if w.exec_():
        pprint.pprint(w.getPointInfo())
    else:
        print("Rejected")

if __name__ == "__main__":
    main()
