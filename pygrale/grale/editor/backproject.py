from PyQt5 import QtCore, QtGui
import numpy as np
from grale.constants import ANGLE_ARCSEC
import openglhelper

def qImageToArray(img):
    img = img.convertToFormat(QtGui.QImage.Format_ARGB32);
    ptr = img.constBits()
    ptr.setsize(img.byteCount())
    arr = np.array(ptr).reshape(img.height(), img.bytesPerLine())
    if img.bytesPerLine() != img.width()*4:
        arr = arr[:,:img.width()*4]
    arr = arr.reshape(img.height(), img.width(), 4)
    arr = arr.astype(np.float)
    return arr

def arrayToQImage(img):
    img = img.copy()
    img = img.astype(np.uint8)
    destImg = QtGui.QImage(img, img.shape[1], img.shape[0], img.shape[1]*4, QtGui.QImage.Format_ARGB32)
    return destImg

def _retraceCPU(ip, img, center, sizes, subSample):

    import grale.images as images

    cx, cy = center
    cx *= ANGLE_ARCSEC
    cy *= ANGLE_ARCSEC
    w, h = sizes
    w *= ANGLE_ARCSEC
    h *= ANGLE_ARCSEC
    arr = qImageToArray(img)

    totalPlane = None
    for idx in range(4):
        data = arr[:,:,idx].reshape(arr.shape[:2])
        src = images.DiscreteSource(data, w, h, [cx, cy])

        plane = ip.renderImages([src], subSamples=subSample)
        if totalPlane is None:
            totalPlane = np.empty(plane.shape + (4,))
        totalPlane[:,:,idx] = plane

    return arrayToQImage(totalPlane)


def backprojectAndRetrace(imgPlane, img, minXY, maxXY, bpXY,useCPU, subSample, numXYGPU):

    centerX, centerY = [ (minXY[i] + maxXY[i])*0.5 for i in range(2) ]
    widthArcsec, heightArcsec = [ (maxXY[i] - minXY[i]) for i in range(2) ]
    ip = imgPlane
    imgSrc, bl, tr = openglhelper.backProject(ip, img, [centerX, centerY], 
                                      [widthArcsec, heightArcsec], [bpXY, bpXY])

    # Recalculate

    w, h = [(tr[i]-bl[i]) for i in [0,1]]
    cx, cy = [ (tr[i]+bl[i])/2 for i in [0,1]]

    if useCPU:
        img = _retraceCPU(ip, imgSrc, [cx, cy], [w, h], subSample)
    else:
        img = imgSrc.mirrored(False, True)
        dims = [ round(numXYGPU*w/h), numXYGPU ] if w < h else [ numXYGPU, round(numXYGPU*h/w) ]
        img = openglhelper.trace(ip, img, [cx, cy], [w, h], dims, subSample)
        img = img.mirrored(False, True)

    return imgSrc, bl, tr, img

