from PyQt5 import QtWidgets, QtCore
import ui_backprojretracedialog
import numpy as np
import pickle
from grale.constants import ANGLE_ARCSEC
import os
from tools import strToBool

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

        settings = QtCore.QSettings()

        def checkBool(key, default=True):
            x = settings.value(key)
            return default if x is None else strToBool(x)

        def check(key, default, conv = None):
            x = settings.value(key)
            if x is None:
                return default
            return x if conv is None else conv(x)

        self.ui.m_splitCheckBox.setChecked(checkBool("bprelensdlg/splitlayers", True))
        self.ui.m_extraBorderEdit.setValue(check("bprelensdlg/extraborder", 0, lambda x : float(x)))
        self.ui.m_imgPixels.setValue(check("bprelensdlg/imgpixels", 512, lambda x : int(x)))
        self.ui.m_bpPixels.setValue(check("bprelensdlg/bppixels", 512, lambda x : int(x)))
        self.ui.m_relensPixels.setValue(check("bprelensdlg/relenspixels", 512, lambda x : int(x)))
        self.ui.m_resample.setValue(check("bprelensdlg/resample", 1, lambda x : int(x)))
        self.ui.m_traceOnlyImages.setChecked(checkBool("bprelensdlg/separateimages", False))
        self.ui.m_bpFileTemplate.setText(check("bprelensdlg/bpfntmpl", _kDefaultBpFileNameTemplate))
        self.ui.m_bpLayerTemplate.setText(check("bprelensdlg/bplayertmpl", _kDefaultBpLayerNameTemplate))
        self.ui.m_relensFileTemplate.setText(check("bprelensdlg/relensfntmpl", _kDefaultRelensFileNameTemplate))
        self.ui.m_relensLayerTemplate.setText(check("bprelensdlg/relenslayertmpl", _kDefaultRelensLayerNameTemplate))

        for edit, default in [ (self.ui.m_bpFileTemplate, _kDefaultBpFileNameTemplate),
                               (self.ui.m_bpLayerTemplate, _kDefaultBpLayerNameTemplate),
                               (self.ui.m_relensFileTemplate, _kDefaultRelensFileNameTemplate),
                               (self.ui.m_relensLayerTemplate, _kDefaultRelensLayerNameTemplate) ]:
            if not edit.text().strip():
                edit.setText(default)
    
    def done(self, r):
        if r == 1 and not self.imagePlane:
            QtWidgets.QMessageBox.warning(self, "No image plane", "No image plane was loaded")
            return

        if r == 1: # Save settings
            settings = QtCore.QSettings()
            settings.setValue("bprelensdlg/splitlayers", self.ui.m_splitCheckBox.isChecked())
            settings.setValue("bprelensdlg/extraborder", self.ui.m_extraBorderEdit.getValue())
            settings.setValue("bprelensdlg/imgpixels", self.ui.m_imgPixels.value())
            settings.setValue("bprelensdlg/bppixels", self.ui.m_bpPixels.value())
            settings.setValue("bprelensdlg/relenspixels", self.ui.m_relensPixels.value())
            settings.setValue("bprelensdlg/resample", self.ui.m_resample.value())
            settings.setValue("bprelensdlg/separateimages", self.ui.m_traceOnlyImages.isChecked())
            settings.setValue("bprelensdlg/bpfntmpl", self.ui.m_bpFileTemplate.text())
            settings.setValue("bprelensdlg/bplayertmpl", self.ui.m_bpLayerTemplate.text())
            settings.setValue("bprelensdlg/relensfntmpl", self.ui.m_relensFileTemplate.text())
            settings.setValue("bprelensdlg/relenslayertmpl", self.ui.m_relensLayerTemplate.text())

        super(BackprojRetraceDialog, self).done(r)
            
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
        return d if d.strip() else os.getcwd()

    def getBackprojectFileTemplate(self):
        t = self.ui.m_bpFileTemplate.text()
        return t if t.strip() else _kDefaultBpFileNameTemplate

    def getBackprojectLayerTemplate(self):
        t = self.ui.m_bpLayerTemplate.text()
        return t if t.strip() else _kDefaultBpLayerNameTemplate

    def getRelensFileTemplate(self):
        t = self.ui.m_relensFileTemplate.text()
        return t if t.strip() else _kDefaultRelensFileNameTemplate

    def getRelensLayerTemplate(self):
        t = self.ui.m_relensLayerTemplate.text()
        return t if t.strip() else _kDefaultRelensLayerNameTemplate

