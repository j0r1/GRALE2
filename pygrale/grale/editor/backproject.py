from PyQt5 import QtCore, QtGui
import numpy as np
from grale.constants import ANGLE_ARCSEC
import openglhelper
import tools
import os

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


def _old_backprojectAndRetrace(imgPlane, img, minXY, maxXY, bpXY,useCPU, subSample, numXYGPU):

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

def _createRGBLayerForImage(img, fn, layerName, bl, tr, yMirror):
    import imagelayer

    if fn is not None:
        img.save(fn)
        l = imagelayer.RGBImageLayer(fn, layerName)
    else:
        l = imagelayer.RGBImageLayer(img, layerName)

    if yMirror:
        mp = [  [ [0, img.height()], bl], 
                [ [0, 0], [ bl[0], tr[1]] ],
                [ [img.width(), img.height()], [tr[0], bl[1]] ],
                [ [img.width(), 0], tr ] ]
    else:
        mp = [  [ [0, 0], bl], 
                [ [0, img.height()], [ bl[0], tr[1]] ],
                [ [img.width(), 0], [tr[0], bl[1]] ],
                [ [img.width(), img.height()], tr ] ]

    l.matchToPoints(mp, True)
    return l

def getImageRegions(scene, layers, splitLayers, extra, numImgPix): # Uses aspect ratio

    exportTimeDelays = False
    exportGroups = False
    imgDat, usedLayers = tools.layersToImagesData(layers, splitLayers, exportGroups, exportTimeDelays, pointsLeftInfo=None, ignoreRemainingPoints=True)

    borders = [ ]
    for i in range(imgDat.getNumberOfImages()):
        border = None
        if imgDat.getNumberOfImagePoints(i) < 3: # Not enough points for triangulation or hull, use first point
            border = [ imgDat.getImagePointPosition(i, 0) ]
        else:
            for fn in [ imgDat.getBorder, imgDat.getConvexHull ]:
                try:
                    border = fn(i)
                except Exception as e:
                    print("Warning:", e)

        if not border:
            raise Exception("Unable to get a border for image {}".format(i+1))
        borders.append(border)

    import grale.images as images
    borders = [ np.array(images.enlargePolygon(b, extra))/ANGLE_ARCSEC for b in borders ]

    imgs = [ ]
    for idx in range(len(borders)):
        border = borders[idx]
        img = scene.getSceneRegionImage_minMax(border.min(0), border.max(0), [ numImgPix, None], border)
        imgs.append(img)

    return list(zip(borders, imgs)), usedLayers

def backprojectAndRetrace(imgPlane, bordersAndImages, 
                       numBPPix = 1024, # same used in x and y direction, is this ok? perhaps we'd lose information otherwise?
                       numRetracePix = 1024, # Again scaled according to aspect ratio
                       numResample = 1,
                       relensSeparately = True,
                       overWriteFiles = False,
                       bpFileNameTemplate = "img_{srcidx}_backproj.png",
                       bpLayerNameTemplate = "Source shape for image {srcidx}: {fn}",
                       relensFileNameTemplate = "img_{srcidx}_to_{tgtidx}_relensed.png",
                       relensLayerNameTemplate = "Relensed source from image {srcidx} to {tgtidx}: {fn}",
                       newImageDir = None,
                       progressCallback = None):

    ip = imgPlane
    numGPUXY = numRetracePix

    if newImageDir is None:
        newImageDir = os.getcwd()
    if progressCallback is None:
        def dummyCb(msg):
            pass
        progressCallback = dummyCb

    newLayers = [ ]

    if not overWriteFiles: # Check what would be overwritten
        filesToWrite = [ ]
        if relensSeparately:
            for idx in range(len(bordersAndImages)):
                filesToWrite.append(os.path.join(newImageDir, bpFileNameTemplate.format(srcidx=idx+1, tgtidx=0)))

                for tgtidx in range(len(bordersAndImages)):
                    filesToWrite.append(os.path.join(newImageDir, relensFileNameTemplate.format(srcidx=idx+1, tgtidx=tgtidx+1)))
        else:
            for idx in range(len(bordersAndImages)):
                for tmpl in [ bpFileNameTemplate, relensFileNameTemplate ]:
                    filesToWrite.append(os.path.join(newImageDir, tmpl.format(srcidx=idx+1, tgtidx=0)))

        filesToOverwrite = [ fn for fn in filesToWrite if os.path.exists(fn) ]
        if filesToOverwrite:
            raise Exception(f"{len(filesToOverwrite)} files would be overwritten, first are:\n" + "\n".join(filesToOverwrite[:10]))

        if len(set(filesToWrite)) != len(filesToWrite):
            raise Exception("Some output files would overwrite each other")

    def addImg(img, srcidx, tgtidx, bl, tr, isSrc):
        yMirror = True if isSrc else False
        fnTemplate = bpFileNameTemplate if isSrc else relensFileNameTemplate
        layerNameTemplate = bpLayerNameTemplate if isSrc else relensLayerNameTemplate

        fn = None if fnTemplate is None else fnTemplate.format(srcidx=srcidx+1, tgtidx=tgtidx+1)
        newLayers.append(_createRGBLayerForImage(
            img, fn, layerNameTemplate.format(srcidx=srcidx+1, tgtidx=tgtidx+1, fn=fn),
            bl, tr, isSrc))

    srcAreas = []
    for idx in range(len(bordersAndImages)):
        border, img = bordersAndImages[idx]

        if ip:
            progressCallback(f"Creating source shape from image {idx+1}")
            centerX, centerY = 0.5*(border.max(0)+border.min(0))
            widthArcsec, heightArcsec = border.max(0)-border.min(0)
            imgSrc, srcbl, srctr = openglhelper.backProject(ip, img, [centerX, centerY], 
                                              [widthArcsec, heightArcsec], [numBPPix, numBPPix])

            srcbl, srctr = np.array(srcbl), np.array(srctr)
        else:
            imgSrc = img.mirrored(False, True)
            srcbl = border.min(0)
            srctr = border.max(0)

        srcCtr = (srcbl+srctr)*0.5
        srcSize = srctr-srcbl

        srcAreas.append([srcbl, srctr])

        addImg(imgSrc, idx, -1, srcbl, srctr, True)

        if numRetracePix <= 0:
            continue

        def getDims(tgtBl, tgtTr):
            tgtW, tgtH = tgtTr-tgtBl 
            return [ round(numGPUXY*tgtW/tgtH), numGPUXY ] if tgtW < tgtH else [ numGPUXY, round(numGPUXY*tgtH/tgtW) ]

        if relensSeparately:
            for tgtidx in range(len(bordersAndImages)):
                tgtBorder = bordersAndImages[tgtidx][0]
                tgtBl, tgtTr = tgtBorder.min(0), tgtBorder.max(0)

                dims = getDims(tgtBl, tgtTr)
                tmpImg = imgSrc.mirrored(False, True)
                
                progressCallback(f"Re-tracing image {tgtidx+1} based on source shape from image {idx+1}")
                tmpImg = openglhelper.trace(ip, tmpImg, srcCtr, srcSize, dims, numResample, tgtBl, tgtTr)
                imgRelens = tmpImg.mirrored(False, True)

                addImg(imgRelens, idx, tgtidx, tgtBl, tgtTr, False)
        else:

            tgtidx = -1

            ri = ip.getRenderInfo()
            tgtBl = np.array(ri["bottomleft"])/ANGLE_ARCSEC
            tgtTr = np.array(ri["topright"])/ANGLE_ARCSEC

            dims = getDims(tgtBl, tgtTr)
            tmpImg = imgSrc.mirrored(False, True)
                
            progressCallback(f"Re-tracing entire image plane based on source shape from image {idx+1}")
            tmpImg = openglhelper.trace(ip, tmpImg, srcCtr, srcSize, dims, numResample, tgtBl, tgtTr)
            imgRelens = tmpImg.mirrored(False, True)

            addImg(imgRelens, idx, tgtidx, tgtBl, tgtTr, False)

    return newLayers, srcAreas

