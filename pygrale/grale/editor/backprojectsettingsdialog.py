from PyQt5 import QtWidgets, QtCore
import ui_backprojectsettingsdialog
import numpy as np
import pickle
from grale.constants import ANGLE_ARCSEC
import os
from tools import strToBool
import helpdialog

class BackprojectSettingsDialog(QtWidgets.QDialog):
    _lastDir = None
    def __init__(self, parent, currentImagePlane):
        super(BackprojectSettingsDialog, self).__init__(parent)

        self.ui = ui_backprojectsettingsdialog.Ui_BackprojectSettingsDialog()
        self.ui.setupUi(self)

        self.ui.m_loadPlaneButton.clicked.connect(self.onLoadImagePlane)
        self.ui.m_inputHelpButton.clicked.connect(self.onInputHelp)
        self.ui.m_outPicHelpButton.clicked.connect(self.onOutputPictureHelp)

        if currentImagePlane is None:
            self.imagePlane = None
            self.ui.m_imgPlaneInfo.setText("No image plane currently loaded")
            self.ui.m_imgPlaneDesc.setText("")
        else:
            self._setImagePlane(currentImagePlane["imgplane"])
            self.ui.m_imgPlaneDesc.setText(currentImagePlane["description"])

        settings = QtCore.QSettings()

        def checkBool(key, default=True):
            x = settings.value(key)
            return default if x is None else strToBool(x)

        def check(key, default, conv = None):
            x = settings.value(key)
            if x is None:
                return default
            return x if conv is None else conv(x)

        self.ui.m_splitCheckBox.setChecked(checkBool("backprojdlg/splitlayers", True))
        self.ui.m_extraBorderEdit.setValue(check("backprojdlg/extraborder", 0, lambda x : float(x)))
        self.ui.m_imgPixels.setValue(check("backprojdlg/imgpixels", 512, lambda x : int(x)))
        self.ui.m_bpPixels.setValue(check("backprojdlg/bppixels", 512, lambda x : int(x)))

    def done(self, r):
        if r == 1 and not self.imagePlane:
            QtWidgets.QMessageBox.warning(self, "No image plane", "No image plane was loaded")
            return

        if r == 1: # Save settings
            settings = QtCore.QSettings()
            settings.setValue("backprojdlg/splitlayers", self.ui.m_splitCheckBox.isChecked())
            settings.setValue("backprojdlg/extraborder", self.ui.m_extraBorderEdit.getValue())
            settings.setValue("backprojdlg/imgpixels", self.ui.m_imgPixels.value())
            settings.setValue("backprojdlg/bppixels", self.ui.m_bpPixels.value())

        super(BackprojectSettingsDialog, self).done(r)
            
    def _setImagePlane(self, ip):
        ri = ip.getRenderInfo()
        tr, bl = [ np.array(ri[x])/ANGLE_ARCSEC for x in [ "topright", "bottomleft"] ]

        self.imagePlane = ip
        self.ui.m_imgPlaneInfo.setText("{}x{} points, ({},{}) -> ({},{}) arcsec".format(
                ri["xpoints"], ri["ypoints"], bl[0], bl[1], tr[0], tr[1]))

    def getImagePlaneInfo(self):
        if not self.imagePlane:
            return None
        return { "imgplane": self.imagePlane, "description": self.ui.m_imgPlaneDesc.text() }

    def onLoadImagePlane(self):
        fn, fltr = QtWidgets.QFileDialog.getOpenFileName(self, "Select image plane file", filter="Pickled image plane files (*.imgplane *.dat)")
        if not fn:
            return

        try:
            ip = pickle.load(open(fn, "rb"))
            ip.getCriticalLines() # check that it's likely an image plane or multi-imageplane
            self._setImagePlane(ip)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, "Can't load image plane", "Error loading image plane: {}".format(e))
    
    def getSplitLayersFlag(self):
        return self.ui.m_splitCheckBox.isChecked()

    def getExtraBorder(self):
        return self.ui.m_extraBorderEdit.getValue()

    def getInputImagePixels(self):
        return self.ui.m_imgPixels.value()

    def getBackProjPixels(self):
        return self.ui.m_bpPixels.value()

    def showHelpDialog(self, msg):
        dlg = helpdialog.HelpDialog(self, msg)
        dlg.exec_()

    def onInputHelp(self):
        self.showHelpDialog("""
        <p><b>Split layer into images</b></p>
        <p>The visible points layers are used to select which parts of the shown scene should be
        projected back onto the source plane. If unchecked, each points layer is interpreted as
        a separate image region. If checked, a points layer may contain more image regions and
        based on triangulation information it will be attempted to split each points layer into
        multiple images.
        </p>
        <p><b>Extra border</b></p>
        <p>The points in an image specify which region will be projected back onto the source
        plane, but it is possible that this way the boundary of the image is not captured
        as desired. With this option, this boundary can be enlarged somewhat.
        </p>
        <p><b>Number of pixels</b></p>
        <p>The rectangular region covering each of the detected images is subdivided into a
        number of pixels. The number specifies the maximum in either X or Y direction,
        the other value will possibly be smaller, to roughly preserve the aspect ratio.
        </p>
        """)

    def onOutputPictureHelp(self):
        self.showHelpDialog("""
        <p><b>Source plane pixels</b></p>
        <p>Each of the image regions that is detected is projected back onto the source plane.
        For each one, a pixelized version is created (and saved to a file) with this amount
        of pixels in both X and Y directions.
        </p>
        """)

