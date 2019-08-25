from PyQt5 import QtWidgets, QtCore, QtGui
from astropy.io import fits
from astropy import wcs
import numpy as np
import uuid
import copy
import grale.images as images # TODO
from grale.constants import * # TODO

from base import Layer, PointGraphicsItemBase, LayerGraphicsItemBase

class ImageGraphicsItem(LayerGraphicsItemBase):
    def __init__(self, layer, parent, drawCrosses):
        super(ImageGraphicsItem, self).__init__(layer, LayerGraphicsItemBase.Cross if drawCrosses else LayerGraphicsItemBase.Circle, True, parent)

        self.imgItem = None

        self.updatePixmap()
        self.updateTransform()

    def updatePixmap(self):
        if not self.imgItem:
            self.imgItem = QtWidgets.QGraphicsPixmapItem(self)

        self.imgItem.setPixmap(self.getLayer().getPixmap())

    def updateTransform(self):
        self.imgItem.setTransform(self.getLayer().getImageTransform())
        self.updatePoints()

class FITSGraphicsItem(ImageGraphicsItem):
    def __init__(self, layer, parent = None):
        super(FITSGraphicsItem, self).__init__(layer, parent, True)

class RGBGraphicsItem(ImageGraphicsItem):
    def __init__(self, layer, parent = None):
        super(RGBGraphicsItem, self).__init__(layer, parent, False)
    
class FITSImageLayer(Layer):
    def __init__(self, fileName, hduIdx = -1, name = "FITS Image layer"):
        super(FITSImageLayer, self).__init__(name)

        f = fits.open(fileName)
        if hduIdx >= 0:
            f0 = f[hduIdx]
        else:
            hduIdx = 0
            while hduIdx < len(f):
                try:
                    f0 = f[hduIdx]
                    s = f0.data.shape
                    if len(s) == 2:
                        break
                except:
                    hduIdx += 1

        #print("Shape =", f0.data.shape)
        self.wcs = wcs.WCS(f0.header)
        self.data = f0.data
        self.fileName = fileName
        self.hduIdx = hduIdx

        self.resetTransform()
        self.resetMinMax()

    def getPixmap(self):
        z0, z1 = min(self.min, self.max), max(self.min, self.max)
        invert = True if self.min > self.max else False
        d = (np.clip(self.data, z0, z1) - z0)/(z1-z0)
        if invert:
            d = 1.0-d
        d = (d * 255.0).astype(np.uint8)
        img = QtGui.QImage(d, self.data.shape[1], self.data.shape[0], self.data.shape[1], QtGui.QImage.Format_Grayscale8)
        return QtGui.QPixmap.fromImage(img)

    def resetMinMax(self):
        self.min = 0
        self.max = 1

    def setMinMax(self, minValue = None, maxValue = None):
        if minValue is not None:
            self.min = minValue
        if maxValue is not None:
            self.max = maxValue

    def getMinMax(self):
        return self.min, self.max

    def resetTransform(self):
        y2, x2 = np.array(self.data.shape)/2.0
        w = self.wcs
        ra, dec = w.all_pix2world(x2, y2, 0)
        self.setCenter(ra, dec)

    def getCenter(self):
        return self.centerRa, self.centerDec
        
    def setCenter(self, ra, dec):
        self.centerRa, self.centerDec = ra, dec
        w = self.wcs

        ra0, dec0 = w.all_pix2world(0, 0, 0)
        ra1, dec1 = w.all_pix2world(self.data.shape[1], 0, 0)
        ra2, dec2 = w.all_pix2world(0, self.data.shape[0], 0)

        xy0, xy1, xy2 = list(map(lambda xy: images.centerOnPosition([xy[0]*ANGLE_DEGREE, xy[1]*ANGLE_DEGREE], 
                                                                    [ra*ANGLE_DEGREE,dec*ANGLE_DEGREE])/ANGLE_ARCSEC,
                                 [(ra0, dec0), (ra1, dec1), (ra2, dec2)]))

        #print(xy0, xy1, xy2)
        Cx, Cy = xy0
        m11 = (xy1[0]-xy0[0])/float(self.data.shape[1])
        m21 = (xy1[1]-xy0[1])/float(self.data.shape[1])
        m12 = (xy2[0]-xy0[0])/float(self.data.shape[0])
        m22 = (xy2[1]-xy0[1])/float(self.data.shape[0])
        
        self.setImageTransform(QtGui.QTransform(m11,m12,m21,m22,Cx,Cy))

    def getPixelCoordinates(self, x, y):
        return self.wcs.all_pix2world(x, y, 0)

    def createGraphicsItem(self):
        return FITSGraphicsItem(self)

    def toSettings(self):
        pts = self.getPoints()

        ptsInfo = [ ]
        for k in pts:
            p = pts[k]
            ptsInfo.append({ "xy": p["xy"], "label": p["label"] })

        settings = {
            "name": self.getName(),
            "fitsfile": self.fileName,
            "hduidx": self.hduIdx,
            "points": ptsInfo,
            "center": [ float(self.centerRa), float(self.centerDec) ],
            "minmax": [ float(self.min), float(self.max) ]
        }
        return settings

    @staticmethod
    def fromSettings(settings):
        l = FITSImageLayer(settings["fitsfile"], settings["hduidx"], settings["name"])
        for p in settings["points"]:
            l.addPoint(p["xy"], p["label"])

        l.setMinMax(*settings["minmax"])
        l.setCenter(*settings["center"])
        return l

class RGBImageLayer(Layer):
    def __init__(self, fileName, name = "RGB Image layer"):
        super(RGBImageLayer, self).__init__(name)

        img = QtGui.QImage()
        if not img.load(fileName):
            raise Exception("Unable to load '{}'".format(fileName))

        self.fileName = fileName
        self.pixmap = QtGui.QPixmap.fromImage(img)

        pixelSize = 0.05
        self.setImageTransform(QtGui.QTransform(-pixelSize, 0, 0, -pixelSize, pixelSize*img.width()/2.0, pixelSize*img.height()/2.0)) # TODO: what's a sensible initialization?

    def getPixmap(self):
        return self.pixmap

    def createGraphicsItem(self):
        return RGBGraphicsItem(self)

    # pointsToMatch should be an array of [ [x,y], [ X, Y] ]
    # entries. The x,y coords refer to the scene coords for the points
    # in this RGB image, the X, Y points are the corresponding scene
    # coordinates to match
    def matchToPoints(self, pointsToMatch, isPixels = False):
        if not pointsToMatch:
            raise Exception("No points to match found")
        
        if len(pointsToMatch) < 3:
            raise Exception("At least three matching points need to be specified")

        target_x = np.array([ XY[0] for xy, XY in pointsToMatch ], dtype=np.double)
        target_y = np.array([ XY[1] for xy, XY in pointsToMatch ], dtype=np.double)
        src_x = np.array([ xy[0] for xy, XY in pointsToMatch ], dtype=np.double)
        src_y = np.array([ xy[1] for xy, XY in pointsToMatch ], dtype=np.double)

        Txx = np.dot(target_x, src_x)
        Txy = np.dot(target_x, src_y)
        Tyx = np.dot(target_y, src_x)
        Tyy = np.dot(target_y, src_y)
        Tx = np.sum(target_x)
        Ty = np.sum(target_y)
        Sxx = np.dot(src_x, src_x)
        Sxy = np.dot(src_x, src_y)
        Syy = np.dot(src_y, src_y)
        Sx = np.sum(src_x)
        Sy = np.sum(src_y)
        N = float(len(pointsToMatch))

        M = np.matrix([
                [ Sxx, Sxy,   0,   0,  Sx,   0 ],
                [ Sxy, Syy,   0,   0,  Sy,   0 ],
                [   0,   0, Sxx, Sxy,   0,  Sx ],
                [   0,   0, Sxy, Syy,   0,  Sy ],
                [  Sx,  Sy,   0,   0,   N,   0 ],
                [   0,   0,  Sx,  Sy,   0,   N ]])


        try:
            MI = M.getI()
        except np.linalg.LinAlgError:
            raise Exception("Can't calculate required transformation, matrix to invert is singular")

        T = np.matrix([
                [ Txx ],
                [ Txy ],
                [ Tyx ],
                [ Tyy ],
                [  Tx ],
                [  Ty ]
            ])

        transformCoefficients = MI*T
        m11, m12, m21, m22, DX, DY = transformCoefficients.A1

        # Note: used a different convention for index numbering, that's
        #       why m21 and m12 are swapped here
        extraTransform = QtGui.QTransform(m11, m21, m12, m22, DX, DY)

        if not isPixels:
            origTransform = self.getImageTransform()
            newTransform = origTransform*extraTransform
        else:
            newTransform = extraTransform

        self.setImageTransform(newTransform)

    def toSettings(self):
        pts = self.getPoints()

        t = self.getImageTransform()
        ptsInfo = [ ]
        for k in pts:
            p = pts[k]
            ptsInfo.append({ "xy": p["xy"], "label": p["label"] })

        settings = {
            "name": self.getName(),
            "rgbfile": self.fileName,
            "points": ptsInfo,
            "transform": [ t.m11(), t.m12(), t.m13(), t.m21(), t.m22(), t.m23(), t.m31(), t.m32(), t.m33() ]
        }
        return settings

    @staticmethod
    def fromSettings(settings):
        l = RGBImageLayer(settings["rgbfile"], settings["name"])
        for p in settings["points"]:
            l.addPoint(p["xy"], p["label"])

        l.setImageTransform(QtGui.QTransform(*settings["transform"]))
        return l
