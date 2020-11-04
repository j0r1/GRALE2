from PyQt5 import QtWidgets, QtCore
import ui_backprojretracedialog
import numpy as np
import pickle
from grale.constants import ANGLE_ARCSEC
import os
from tools import strToBool
import helpdialog

_kDefaultBpFileNameTemplate = "img_{srcidx}_backproj.png"
_kDefaultBpLayerNameTemplate = "Source shape for image {srcidx}: {fn}"
_kDefaultRelensFileNameTemplate = "img_{srcidx}_to_{tgtidx}_relensed.png"
_kDefaultRelensLayerNameTemplate = "Relensed source from image {srcidx} to {tgtidx}: {fn}"

class BackprojRetraceDialog(QtWidgets.QDialog):
    _lastDir = None
    def __init__(self, parent, currentImagePlane):
        super(BackprojRetraceDialog, self).__init__(parent)

        self.ui = ui_backprojretracedialog.Ui_BackprojRetraceDialog()
        self.ui.setupUi(self)

        self.ui.m_loadPlaneButton.clicked.connect(self.onLoadImagePlane)
        self.ui.m_inputHelpButton.clicked.connect(self.onInputHelp)
        self.ui.m_outPicHelpButton.clicked.connect(self.onOutputPictureHelp)
        self.ui.m_templateHelpButton.clicked.connect(self.onTemplateHelp)
        self.ui.m_browseButton.clicked.connect(self.onBrowseDirectory)

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

        if BackprojRetraceDialog._lastDir:
            self.ui.m_outputDirectory.setText(BackprojRetraceDialog._lastDir)
        if not self.ui.m_outputDirectory.text().strip():
            self.ui.m_outputDirectory.setText(os.getcwd())
    
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
            BackprojRetraceDialog._lastDir = self.ui.m_outputDirectory.text()

        super(BackprojRetraceDialog, self).done(r)
            
    def _setImagePlane(self, ip):
        ri = ip.getRenderInfo()
        tr, bl = [ np.array(ri[x])/ANGLE_ARCSEC for x in [ "topright", "bottomleft"] ]

        self.imagePlane = ip
        self.ui.m_imgPlaneInfo.setText("{}x{} points, ({:.2f},{:.2f}) -> ({:.2f},{:.2f}) arcsec".format(
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
        <p><b>Re-lens grid</b></p>
        <p>When re-lensing the source shapes that are recovered this way, this specifies the
        amount of pixels used for the calculation. This is again the maximum in either X or Y
        direction, the other direction may be smaller, depending on the aspect ratio. The
        precise meaning of the grid depends on the 'trace only...' option below.
        </p>
        <p><b>Subsample</b></p>
        <p>
        If this is 1, the center of each grid pixel is traced to the source plane to get the
        corresponding source color. If more than one, each pixel is subdivided this many times
        in both X and Y directions and the resulting pixel value is the average of the colors
        at the corresponding source plane positions.
        </p>
        <p><b>Trace only the image regions</b></p>
        <p>If unchecked, the entire image plane region, as defined in the image plane file
        that has been loaded, will be recalculated. The resolution of the grid above then
        specifies the resolution of this entire region. The image regions that were used
        as input will be recovered, but probably with some loss of detail as they are smaller
        than the entire image plane. Calculating the entire image plane is of course a good
        way to get an overall view of all the generated images.
        </p>
        <p>If checked, then only the same image regions that are used as input are recalculated,
        each one being covered with a grid of the specified resolution. This will improve the
        detail in each individual image, but of course may miss some features as not the entire
        image plane is recalculated. For each backproject image, all images are recalculated,
        so this will require more iterations as the 'unchecked' version.
        </p>
        """)

    def onTemplateHelp(self):
        self.showHelpDialog("""
        <p><b>Back-projected image: file name</b></p>
        <p>Each back-projected image region is saved to a file, this specified the template
        file name. The '{srcidx}' part will be replaced by the number of the image region.
        </p>
        <p><b>Back-projected image: layer name</b></p>
        <p>For each of the previous files, a new layer will be added to the editor, and this
        specifies the template of this layer's name. Again '{srcidx}' can be used to indicate
        which image region it is based on. Here, '{fn}' can be used to substitute the 
        file name above.
        </p>
        <p><b>Re-lensed image plane/image: file name</b></p>
        <p>Here too, the generated pixelized regions are saved to a file. As before, '{srcidx}' indicates
        which back-projected image is used as the source. If the entire image plane is calculated
        at once, then '{tgtidx}' is 0, but if each image region is calculated separately, this
        indicates which region is calculated.
        </p>
        <p><b>Re-lensed image plane/image: layer name</b></p>
        <p>Similar as before, a layer will be added for each of these files, and '{srcidx}' and
        '{tgtidx}' have the same meaning. The identifier '{fn}' can again be used to get the
        corresponding file name.
        </p>
        """)

    def onBrowseDirectory(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(self, "Select directory to store new image files in",
                                                       self.getOutputDirectory(), QtWidgets.QFileDialog.ShowDirsOnly|QtWidgets.QFileDialog.DontResolveSymlinks)
        if not d:
            return
        self.ui.m_outputDirectory.setText(d)

