from PyQt5 import QtWidgets, QtCore
import ui_imgregionsettingsdialog
import numpy as np
import pickle
from grale.constants import ANGLE_ARCSEC
import os
from tools import strToBool
import helpdialog

class ImageRegionSettingsDialog(QtWidgets.QDialog):
    def __init__(self, parent):
        super(ImageRegionSettingsDialog, self).__init__(parent)

        self.ui = ui_imgregionsettingsdialog.Ui_ImageRegionSettingsDialog()
        self.ui.setupUi(self)

        self.ui.m_inputHelpButton.clicked.connect(self.onInputHelp)

        settings = QtCore.QSettings()

        def checkBool(key, default=True):
            x = settings.value(key)
            return default if x is None else strToBool(x)

        def check(key, default, conv = None):
            x = settings.value(key)
            if x is None:
                return default
            return x if conv is None else conv(x)

        self.ui.m_splitCheckBox.setChecked(checkBool("imgregdlg/splitlayers", True))
        self.ui.m_extraBorderEdit.setValue(check("imgregdlg/extraborder", 0, lambda x : float(x)))
        self.ui.m_imgPixels.setValue(check("imgregdlg/imgpixels", 512, lambda x : int(x)))

    def done(self, r):
        if r == 1: # Save settings
            settings = QtCore.QSettings()
            settings.setValue("imgregdlg/splitlayers", self.ui.m_splitCheckBox.isChecked())
            settings.setValue("imgregdlg/extraborder", self.ui.m_extraBorderEdit.getValue())
            settings.setValue("imgregdlg/imgpixels", self.ui.m_imgPixels.value())

        super(ImageRegionSettingsDialog, self).done(r)
            
    def getSplitLayersFlag(self):
        return self.ui.m_splitCheckBox.isChecked()

    def getExtraBorder(self):
        return self.ui.m_extraBorderEdit.getValue()

    def getInputImagePixels(self):
        return self.ui.m_imgPixels.value()

    def showHelpDialog(self, msg):
        dlg = helpdialog.HelpDialog(self, msg)
        dlg.exec_()

    def onInputHelp(self):
        self.showHelpDialog("""
        <p><b>Split layer into images</b></p>
        <p>The visible points layers are used to select which parts of the shown scene should be
        shown as image regions. If unchecked, each points layer is interpreted as
        a separate image region. If checked, a points layer may contain more image regions and
        based on triangulation information it will be attempted to split each points layer into
        multiple images.
        </p>
        <p><b>Extra border</b></p>
        <p>The points in an image specify the region that will be shown, but it is
        possible that this way the boundary of the image is not captured
        as desired. With this option, this boundary can be enlarged somewhat.
        </p>
        <p><b>Number of pixels</b></p>
        <p>The rectangular region covering each of the detected images is subdivided into a
        number of pixels. The number specifies the maximum in either X or Y direction,
        the other value will possibly be smaller, to roughly preserve the aspect ratio.
        </p>
        """)

