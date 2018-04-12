from PyQt5 import QtWidgets, QtCore, QtGui
from astropy.io import fits
from astropy import wcs
import numpy as np
import uuid
import copy
import grale.images as images # TODO
from grale.constants import * # TODO

from base import Layer, EmptyGraphicsItem, PointGraphicsItemBase

class ImageLayerBase(Layer):
    def __init__(self, name):
        super(ImageLayerBase, self).__init__(name)

    def _setPointKw(self, pointDict, label):
        pointDict["label"] = label

    def getPixmap(self):
        raise Exception("TODO: ImageLayerBase.getPixmap needs to be implemented in derived class")

class MatchPointGraphicsItem(PointGraphicsItemBase):

    s_normalLinePen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.red), 0.2)
    s_selectedLinePen = QtGui.QPen(QtGui.QBrush(QtCore.Qt.blue), 0.2)

    s_rectBrush = QtGui.QBrush(QtGui.QColor(0, 0, 0, 1))
    s_rectPen = QtGui.QPen(s_rectBrush, 0.01)

    def __init__(self, parent, uuid, label, cross):
        super(MatchPointGraphicsItem, self).__init__(parent, uuid, label)

        self.parent = parent

        rectBrush = MatchPointGraphicsItem.s_rectBrush
        rectPen = MatchPointGraphicsItem.s_rectPen
        self.rect = QtWidgets.QGraphicsRectItem(-1.1, -1.1, 2.2, 2.2, self)
        self.rect.setPen(rectPen)
        self.rect.setBrush(rectBrush)

        linePen = MatchPointGraphicsItem.s_normalLinePen
        if cross:
            self.line1 = QtWidgets.QGraphicsLineItem(-1,-1, 1, 1, self)
            self.line1.setPen(linePen)

            self.line2 = QtWidgets.QGraphicsLineItem(-1, 1, 1,-1, self)
            self.line2.setPen(linePen)

            self.selectionChangeItems = [ self.line1, self.line2 ]
        else:
            self.circle = QtWidgets.QGraphicsEllipseItem(-1, -1, 2, 2, self)
            self.circle.setPen(linePen)

            self.selectionChangeItems = [ self.circle ]

    def scenePositionToPointPosition(self, pos):
        return self.parent.getPixelPosition(pos)

    def pointPositionToScenePosition(self, pos):
        return self.parent.getScenePosition(pos)

    def onSelected(self, value):
        p = MatchPointGraphicsItem.s_selectedLinePen if value else MatchPointGraphicsItem.s_normalLinePen
        for i in self.selectionChangeItems:
            i.setPen(p)

class ImageGraphicsItem(EmptyGraphicsItem):
    def __init__(self, layer, parent, drawCrosses):
        super(ImageGraphicsItem, self).__init__(parent)

        self.drawCrosses = drawCrosses
        self.layer = layer
        self.imgItem = None
        self.matchPointItems = [ ]

        self.updatePixmap()
        self.updateTransform()

    def arePointsVisible(self):
        s = self.scene()
        if not s:
            return True
        return s.areMatchPointsVisible()

    def getLayer(self):
        return self.layer

    def getScenePosition(self, pos):
        p = self.imgItem.mapToScene(pos[0], pos[1])
        return (p.x(), p.y())

    def getPixelPosition(self, pos):
        p = self.imgItem.mapFromScene(pos[0], pos[1])
        return (p.x(), p.y())

    def addPoint(self, pos, label = None, uuid = None):
        pixPos = self.getPixelPosition(pos)
        if not uuid:
            uuid = self.layer.addPoint(pixPos, label=label)
        else:
            self.layer.setPoint(uuid, pixPos, label=label)
        self.updatePoints()
        return uuid

    def _addMatchPointItem(self, uuid, label, pos):
        # TODO: check if uuid already present as safety check?

        m = MatchPointGraphicsItem(self, uuid, label, self.drawCrosses)
        m.setPos(pos[0], pos[1])
        m.setVisible(self.arePointsVisible())

        self.matchPointItems.append(m)

    def updatePoints(self):
        
        # clear existing points
        s = self.scene()
        for i in self.matchPointItems:
            if s: s.removeItem(i)

        self.matchPointItems = [ ]

        pts = self.layer.getPoints()
        for k in pts:
            pixelPos = pts[k]["xy"]
            p = self.imgItem.mapToScene(pixelPos[0], pixelPos[1])
            self._addMatchPointItem(k, pts[k]["label"], (p.x(), p.y()))

    def updatePixmap(self):
        if not self.imgItem:
            self.imgItem = QtWidgets.QGraphicsPixmapItem(self)

        self.imgItem.setPixmap(self.layer.getPixmap())

    def updateTransform(self):
        self.imgItem.setTransform(self.layer.getTransform())
        self.updatePoints()

class FITSGraphicsItem(ImageGraphicsItem):
    def __init__(self, layer, parent = None):
        super(FITSGraphicsItem, self).__init__(layer, parent, True)

class RGBGraphicsItem(ImageGraphicsItem):
    def __init__(self, layer, parent = None):
        super(RGBGraphicsItem, self).__init__(layer, parent, False)
    
class FITSImageLayer(ImageLayerBase):
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
        d = ((np.clip(self.data, self.min, self.max) - self.min)/(self.max-self.min) * 255.0).astype(np.uint8)
        img = QtGui.QImage(d, self.data.shape[1], self.data.shape[0], QtGui.QImage.Format_Grayscale8)
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
        
        self.setTransform(QtGui.QTransform(m11,m12,m21,m22,Cx,Cy))

    def getPixelCoordinates(self, x, y):
        return self.wcs.all_pix2world(x, y, 0)

    def createGraphicsItem(self):
        return FITSGraphicsItem(self)

    def toSettings(self):
        pts = self.getPoints()

        settings = {
            "name": self.getName(),
            "fitsfile": self.fileName,
            "hduidx": self.hduIdx,
            "points": [ pts[k] for k in pts ],
            "center": [ float(self.centerRa), float(self.centerDec) ],
            "minmax": [ float(self.min), float(self.max) ]
        }
        return settings

    @staticmethod
    def fromSettings(settings):
        l = FITSImageLayer(settings["fitsfile"], settings["hduidx"], settings["name"])
        for p in settings["points"]:
            p = copy.deepcopy(p)
            xy = p["xy"]
            del p["xy"]
            l.addPoint(xy, **p)

        l.setMinMax(*settings["minmax"])
        l.setCenter(*settings["center"])
        return l

class RGBImageLayer(ImageLayerBase):
    def __init__(self, fileName, name = "RGB Image layer"):
        super(RGBImageLayer, self).__init__(name)

        img = QtGui.QImage()
        if not img.load(fileName):
            raise Exception("Unable to load '{}'".format(fileName))

        self.fileName = fileName
        self.pixmap = QtGui.QPixmap.fromImage(img)

        pixelSize = 0.05
        self.setTransform(QtGui.QTransform(-pixelSize, 0, 0, -pixelSize, pixelSize*img.width()/2.0, pixelSize*img.height()/2.0)) # TODO: what's a sensible initialization?

    def getPixmap(self):
        return self.pixmap

    def createGraphicsItem(self):
        return RGBGraphicsItem(self)

    # pointsToMatch should be an array of [ [x,y], [ X, Y] ]
    # entries. The x,y coords refer to the scene coords for the points
    # in this RGB image, the X, Y points are the corresponding scene
    # coordinates to match
    def matchToPoints(self, pointsToMatch):
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
        origTransform = self.getTransform()
        newTransform = origTransform*extraTransform

        self.setTransform(newTransform)

    def toSettings(self):
        pts = self.getPoints()

        t = self.getTransform()
        settings = {
            "name": self.getName(),
            "rgbfile": self.fileName,
            "points": [ pts[k] for k in pts ],
            "transform": [ t.m11(), t.m12(), t.m13(), t.m21(), t.m22(), t.m23(), t.m31(), t.m32(), t.m33() ]
        }
        return settings

    @staticmethod
    def fromSettings(settings):
        l = RGBImageLayer(settings["rgbfile"], settings["name"])
        for p in settings["points"]:
            p = copy.deepcopy(p)
            xy = p["xy"]
            del p["xy"]
            l.addPoint(xy, **p)

        l.setTransform(QtGui.QTransform(*settings["transform"]))
        return l
