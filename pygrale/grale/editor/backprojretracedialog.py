from PyQt5 import QtWidgets, QtCore
import ui_backprojretracedialog
import numpy as np
import pickle
from grale.constants import ANGLE_ARCSEC
import os

_kDefaultBpFileNameTemplate = "img_{srcidx}_backproj.png"
_kDefaultBpLayerNameTemplate = "Source shape for image {srcidx}: {fn}"
_kDefaultRelensFileNameTemplate = "img_{srcidx}_to_{tgtidx}_relensed.png"
_kDefaultRelensLayerNameTemplate = "Relensed source from image {srcidx} to {tgtidx}: {fn}"

class BackprojRetraceDialog(QtWidgets.QDialog):
    def __init__(self, parent, currentImagePlane):
        super(BackprojRetraceDialog, self).__init__(parent)

        self.ui = ui_backprojretracedialog.Ui_BackprojRetraceDialog()
        self.ui.setupUi(self)

        self.ui.m_loadPlaneButton.clicked.connect(self.onLoadImagePlane)

        if currentImagePlane is None:
            self.imagePlane = None
            self.ui.m_imgPlaneInfo.setText("No image plane currently loaded")
            self.ui.m_imgPlaneDesc.setText("")
        else:
            self._setImagePlane(currentImagePlane["imgplane"])


            self.ui.m_imgPlaneDesc.setText(currentImagePlane["description"])
            
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

    def getRelensPixels(self):
        return self.ui.m_relensPixels.value()

    def getResampleParameter(self):
        return self.ui.m_resample.value()

    def getLensSeparateImages(self):
        return self.ui.m_traceOnlyImages.isChecked()

    def getOverwriteFlag(self):
        return self.ui.m_overwriteCheckBox.isChecked()

    def getOutputDirectory(self):
        d = self.ui.m_outputDirectory.text()
        return d if d else os.getcwd()

    def getBackprojectFileTemplate(self):
        t = self.ui.m_bpFileTemplate.text()
        return t if t else _kDefaultBpFileNameTemplate

    def getBackprojectLayerTemplate(self):
        t = self.ui.m_bpLayerTemplate.text()
        return t if t else _kDefaultBpLayerNameTemplate

    def getRelensFileTemplate(self):
        t = self.ui.m_relensFileTemplate.text()
        return t if t else _kDefaultRelensFileNameTemplate

    def getRelensLayerTemplate(self):
        t = self.ui.m_relensLayerTemplate.text()
        return t if t else _kDefaultRelensLayerNameTemplate

