from PyQt5 import QtWidgets
import ui_fitslayerinfodialog

class FITSLayerInfoDialog(QtWidgets.QDialog):
    def __init__(self, centerRaDec, minMax, parent = None):
        super(FITSLayerInfoDialog, self).__init__(parent)

        self.ui = ui_fitslayerinfodialog.Ui_FITSLayerInfoDialog()
        self.ui.setupUi(self)

        self.ui.m_raEdit.setValue(centerRaDec[0])
        self.ui.m_decEdit.setValue(centerRaDec[1])
        self.ui.m_minEdit.setValue(minMax[0])
        self.ui.m_maxEdit.setValue(minMax[1])

    def getCenter(self):
        return ( self.ui.m_raEdit.getValue(), self.ui.m_decEdit.getValue() )

    def getMinMax(self):
        return ( self.ui.m_minEdit.getValue(), self.ui.m_maxEdit.getValue() )

