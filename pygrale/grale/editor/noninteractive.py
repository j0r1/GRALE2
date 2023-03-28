from base import GraphicsView, Layer
import scenes
import sys
import json
import pprint
import multiprocessing
import checkqt
import tools
import multiprocessing
import copy
import os
import backproject
import numpy as np
from grale.constants import ANGLE_ARCSEC

class SceneViewDesciptionException(Exception):
    pass

class SceneViewDesciption(object):
    def __init__(self, json_obj_fn = None):
        if json_obj_fn is None:
            self.clear()
            return

        if type(json_obj_fn) == dict:
            self.loadObject(json_obj_fn)
        elif os.path.exists(json_obj_fn):
            self.loadFile(json_obj_fn)
        else:
            self.loadJson(json_obj_fn)

    def clear(self):
        self.activeLayerIndex = -1

        self.globalSettings = {
            "axisleft": False,
            "ispointsizefixed": True,
            "nextmatchpointfits": 0,
            "nextmatchpointrgb": 0,
            "pointsizearcsecsize": 1,
            "pointsizepixelsize": 16,
            "showaxes": True,
            "showmatchpoints": True
        }
        
        self.layers = []
        self.visibilities = []

        self.view = { "center": [ 0.0, 0.0 ], "zoom": 1.0 }

    def getLayerType(self, i):
        if "rgbfile" in self.layers[i]:
            return "rgb"
        if "fitsfile" in self.layers[i]:
            return "fits"
        return "points"

    def getNumberOfLayers(self):
        return len(self.layers)

    def getLayerName(self, i):
        return self.layers[i]["name"]

    def isLayerVisible(self, i):
        return self.visibilities[i]

    def setLayerVisibility(self, i, v):
        self.visibilities[i] = bool(v)

    def setAxisVisible(self, v):
        self.globalSettings["showaxes"] = bool(v)

    def isAxisVisible(self):
        return self.globalSettings["showaxes"]

    def loadFile(self, fn):
        obj = json.load(open(fn, "rt"))
        self.loadObject(obj)

    def loadJson(self, jsonText):
        obj = json.loads(jsonText)
        self.loadObject(obj)

    def loadObject(self, obj):
        self.clear()

        if "view" in obj:
            self.view["zoom"] = float(obj["view"]["zoom"])
            self.view["center"] = [ float(obj["view"]["center"][0]), float(obj["view"]["center"][1]) ]
        
        layersToAdd = None
        visibilitiesToAdd = None

        if "state" in obj:
            layersToAdd = obj["state"]["layers"] # the 'addLayers' function will make the copy
            visibilitiesToAdd = obj["state"]["visibilities"]

            if len(visibilitiesToAdd) != len(layersToAdd):
                raise SceneViewDesciptionException("Number of visibilities doesn't match number of layers")

            self.activeLayerIndex = int(obj["state"]["activelayerindex"])

            gbs = obj["state"]["globalsettings"]
            self.globalSettings["axisleft"] = bool(gbs["axisleft"])
            self.globalSettings["ispointsizefixed"] = bool(gbs["ispointsizefixed"])
            self.globalSettings["nextmatchpointfits"] = int(gbs["nextmatchpointfits"])
            self.globalSettings["nextmatchpointrgb"] = int(gbs["nextmatchpointrgb"])
            self.globalSettings["pointsizearcsecsize"] = float(gbs["pointsizearcsecsize"])
            self.globalSettings["pointsizepixelsize"] = int(gbs["pointsizepixelsize"])
            self.globalSettings["showaxes"] = bool(gbs["showaxes"])
            self.globalSettings["showmatchpoints"] = bool(gbs["showmatchpoints"])

        else: # Assume it's layers only
            layersToAdd = [ l for i in obj ] # the 'addLayers' function will make the copy
            visibilitiesToAdd = [ True for _ in layersToAdd ] # Set everything to visible

        self.addLayers(layersToAdd, visibilitiesToAdd)

    def addLayers(self, layers, visibilities = None):
        if visibilities is None:
            visibilities = [ True for _ in layers ]

        if len(layers) != len(visibilities):
            raise SceneViewDesciptionException("Not the same number of layers and visibilities")

        for l, v in zip(layers, visibilities):
            self.addLayer(l, v)

    def addLayer(self, layer, visibility = True):
        # Note: if we run this as check, we need to have at least a QGuiApplication, and the
        #       idea is to only need that in a background process
        # Layer.fromSettings(layer) # For now we can't

        # Do other elementary checks:
        checkTd = False
        if "rgbfile" in layer:
            if not "transform" in layer:
                raise SceneViewDesciptionException("Expecting 'transform' in RGB layer")
        elif "fitsfile" in layer:
            for k in [ "hduidx", "center", "minmax" ]:
                if not k in layer:
                    raise SceneViewDesciptionException("Expecting '{}' in FITS layer".format(k))
        else:
            if not "triangles" in layer:
                raise SceneViewDesciptionException("Expecting a 'triangles' entry in points layer")
            checkTd = True

        # Check name
        if not "name" in layer:
            raise Exception("Expecting a 'name' for the layer")

        for pt in layer["points"][:1]: # TODO: to save time, we'll only check the first point?
            if not "label" in pt:
                raise Exception("Expecting a 'label' for each point")
            if not "xy" in pt:
                raise Exception("Exception a 'xy' coordinate for each point")
            if checkTd and not "timedelay" in pt:
                raise Exception("Expecting a 'timedelay' entry for each point")

        self.layers.append(copy.deepcopy(layer))
        self.visibilities.append(bool(visibility))

    def toObject(self):
        obj = {
            "state": {
                "activelayerindex": self.activeLayerIndex,
                "globalsettings": copy.deepcopy(self.globalSettings),
                "layers": copy.deepcopy(self.layers),
                "visibilities": copy.deepcopy(self.visibilities)
            },
            "view": copy.deepcopy(self.view)
        }
        return obj

    def getLayer(self, i):
        return copy.deepcopy(self.layers[i])

    def toLayerObject(self, i):
        return Layer.fromSettings(self.layers[i])

    def saveFile(self, fn):
        json.dump(self.toObject(), open(fn, "wt"), indent=2)

class MultipleLayerScene(scenes.LayerScene):

    def __init__(self, sceneViewDesc = None):
        super(MultipleLayerScene, self).__init__()

        self.layersAndWidgets = []
        self.nextZ = 0

        if sceneViewDesc:
            self.loadSceneViewDescription(sceneViewDesc)

    def addLayer(self, l):
        w = l.createGraphicsItem()
        self.layersAndWidgets.append((l, w))

        self.nextZ += 1
        w.setZValue(self.nextZ)
        self.addItem(w)
        w.updatePoints()

#    def getLayerItem(self, uuid):
#        if uuid in self.layerItems:
#            return self.layerItems[uuid]
#        raise Exception("Specified layer {} not found".format(uuid))

    def getCurrentItemAndLayer(self):
        return self.layersAndWidgets[-1]

    def loadSceneViewDescription(self, svd):
        d = svd.globalSettings
        self.setAxesVisible(d["showaxes"])
        self.setAxisLeft(d["axisleft"])
        self.setPointSizePixelSize(d["pointsizepixelsize"])
        self.setPointSizeSceneSize(d["pointsizearcsecsize"])
        self.setPointSizeFixed(d["ispointsizefixed"])
        self.setMatchPointsVisible(d["showmatchpoints"])

        scale = svd.view["zoom"]
        self.onScaleChanged(scale) # May be needed to recalculate the point size transformation

        for layer,visible in zip(svd.layers, svd.visibilities):
            if visible:
                self.addLayer(Layer.fromSettings(layer))

def _background_helper(q1, q2, func):
    params = q1.get()

    try:
        from PyQt5 import QtWidgets, QtCore, QtGui
        app = QtWidgets.QApplication(sys.argv)

        res = func(params)
        q2.put([True, res])
    except Exception as e:
        q2.put([False, str(e)])

def _background_helper_getSceneRegionImageNumPy(q1, q2):
    _background_helper(q1, q2, lambda params: MultipleLayerScene(SceneViewDesciption(params[0])).getSceneRegionNumPyArray(*params[1:]))

def _foreground_helper(helperFunction, paramList):

    ctx = None
    # The default 'fork' on Unix doesn't work nice with Qt
    try:
        ctx = multiprocessing.get_context("forkserver") # TODO: is this better if available?
    except Exception as e:
        pass

    if not ctx:
        ctx = multiprocessing.get_context("spawn")

    checkqt.checkQtAvailable()

    q1, q2 = ctx.Queue(), ctx.Queue()
    p = ctx.Process(target=helperFunction, args=(q1,q2))
    p.start()
    q1.put(paramList)
    flag, result = q2.get()
    p.join()
    if not flag:
        raise Exception(result)
    return result

def getSceneRegionNumPyArray(sceneViewDesc, bottomLeft, topRight, widthPixels = None, heightPixels = None, grayScale = False):
    return _foreground_helper(_background_helper_getSceneRegionImageNumPy,
                              [ sceneViewDesc.toObject(), 
                                [ bottomLeft[0]/ANGLE_ARCSEC, bottomLeft[1]/ANGLE_ARCSEC], 
                                [ topRight[0]/ANGLE_ARCSEC, topRight[1]/ANGLE_ARCSEC], 
                                widthPixels, heightPixels, grayScale ])

def _background_helper_getPointsLayersFromImagesData(q1, q2):
    _background_helper(q1, q2, lambda params: [ l.toSettings() for l in tools.importImagesDataToLayers(*params) ])

def getPointsLayersFromImagesData(imgDat, which, layerTitleFormat = "Points from images data"):
    return _foreground_helper(_background_helper_getPointsLayersFromImagesData, [ imgDat, which, layerTitleFormat ])

def _background_helper_getImagesDataFromVisibleLayers(q1, q2):
    _background_helper(q1, q2, lambda params: tools.layersToImagesData([ Layer.fromSettings(l) for l in params[0] ], *params[1:])[0])

def getImagesDataFromVisibleLayers(sceneViewDesc, splitLayers = True, exportGroups = True, exportTimeDelays = True):
    layers = [ l for l,v in zip(sceneViewDesc.layers, sceneViewDesc.visibilities) if v ]
    return _foreground_helper(_background_helper_getImagesDataFromVisibleLayers, [ layers, splitLayers, exportGroups, exportTimeDelays ])

def _background_helper_getLayersFromBackProjectRetrace(q1, q2):
    def f(params):
        svdObj = params[0]
        imagePlane = params[1]
        splitLayers = params[2]
        extra = params[3]
        numPix = params[4]

        visiblePointLayers = []

        # Hide points layers to get image regions
        svd = SceneViewDesciption(svdObj)
        origVis = [ svd.isLayerVisible(i) for i in range(svd.getNumberOfLayers()) ]
        for i in range(svd.getNumberOfLayers()):
            if svd.getLayerType(i) == "points":
                if svd.isLayerVisible(i):
                    visiblePointLayers.append(svd.toLayerObject(i))
                    svd.setLayerVisibility(i, False)

        sceneNoPts = MultipleLayerScene(svd)
        bordersAndImages, usedLayers = backproject.getImageRegions(sceneNoPts, visiblePointLayers, splitLayers, extra, numPix)    

        # Restore points layers
        for i, v in enumerate(origVis):
            svd.setLayerVisibility(i, v)

        #def cb(s):
        #    print(s)

        newLayers, srcAreas = backproject.backprojectAndRetrace(imagePlane, bordersAndImages, *params[5:]) #, cb)
        newLayers = [ l.toSettings() for l in newLayers ]
        srcAreas = [ [ c.astype(np.float64)*ANGLE_ARCSEC for c in s ] for s in srcAreas ]
        return newLayers, srcAreas

    _background_helper(q1, q2, f)

# TODO
#  - sceneViewDesc
#  - ip: image plane instance
#  - splitLayers: split layer into multiple images
#  - extra: extra border to add
#  - numPix: image regions are covered by this number of pixels, aspect ratio is used
#  - numBPPix: backprojected image is used as source, with numBPPix*numBPPix pixels
#  - numRetracePix: if relensSeparately is False, the entire image plane region will be covered by, set to negative to backproject only
#  - numResample: ?
#  - relensSeparately: recalculate image regions only, based on other backprojected images
#  - overWriteFiles:
#  - bpFileNameTemplate:
#  - bpLayerNameTemplate:
#  - relensFileNameTemplate:
#  - relensLayerNameTemplate:
#  - newImageDir
def getLayersFromBackProjectRetrace(sceneViewDesc, ip,
        splitLayers = True,
        extra = 0.1*ANGLE_ARCSEC,
        numPix = 1024,
        numBPPix = 1024, # same used in x and y direction, is this ok? perhaps we'd lose information otherwise?
        numRetracePix = 1024, # Again scaled according to aspect ratio , set to negative to backproject only!
        numResample = 1,
        relensSeparately = True,
        overWriteFiles = False,
        bpFileNameTemplate = "img_{srcidx}_backproj.png",
        bpLayerNameTemplate = "Source shape for image {srcidx}: {fn}",
        relensFileNameTemplate = "img_{srcidx}_to_{tgtidx}_relensed.png",
        relensLayerNameTemplate = "Relensed source from image {srcidx} to {tgtidx}: {fn}",
        newImageDir = None):

    return _foreground_helper(_background_helper_getLayersFromBackProjectRetrace,
            [ sceneViewDesc.toObject(), 
              ip, splitLayers, extra, numPix, numBPPix, numRetracePix, numResample, relensSeparately,
                                      overWriteFiles, bpFileNameTemplate, bpLayerNameTemplate, relensFileNameTemplate,
                                      relensLayerNameTemplate, newImageDir ])

def getSourceShapeFromSourceRGBLayer(layer, channel="gray"): # or red, green, blue
    transform = layer["transform"]

    import matplotlib.pyplot as plt
    img = plt.imread(layer["rgbfile"])
    assert len(img.shape) == 3 and img.shape[2] >= 3
    
    xPixScale = transform[0]
    yPixScale = transform[4]
    xOffsetArcsec = transform[6]
    yOffsetArcsec = transform[7]

    assert abs(transform[1]) < 1e-6 and abs(transform[2]) < 1e-6 and abs(transform[3]) < 1e-6 and abs(transform[5]) < 1e-6 and abs(transform[8]) == 1, "Doesn't seem to be a transformation for a source shape"

    heightPixels, widthPixels = img.shape[:2]
    heightArcsec = yPixScale*heightPixels
    widthArcsec = xPixScale*widthPixels

    print(widthArcsec, heightArcsec)

    assert widthArcsec > 0 and heightArcsec < 0

    xCtrArcsec = xOffsetArcsec + widthArcsec/2.0
    yCtrArcsec = yOffsetArcsec + heightArcsec/2.0
    ctr = np.array([xCtrArcsec, yCtrArcsec])*ANGLE_ARCSEC

    if channel == "gray":
        imgPartR, imgPartG, imgPartB = img[:,:,0], img[:,:,1], img[:,:,2]
        values = 0.299*imgPartR + 0.587*imgPartG + 0.114*imgPartB
    elif channel == "red":
        values = img[:,:,0]
    elif channel == "green":
        values = img[:,:,1]
    elif channel == "blue":
        values = img[:,:,2]
    else:
        raise Exception("Invalid 'channel' setting")

    import grale.images as images
    src = images.DiscreteSource(values.astype(np.float64), abs(widthArcsec)*ANGLE_ARCSEC, abs(heightArcsec)*ANGLE_ARCSEC, ctr)

    return src

def main():

    #svd = SceneViewDesciption("cl0024.json")
    #pprint.pprint(svd.toObject())

    #svd = SceneViewDesciption("cl0024.json")
    #img = getSceneRegionNumPyArray(svd, [ -10, -10 ], [10, 10], 512)
    #
    #import matplotlib.pyplot as plt
    #plt.imshow(img)
    #plt.show()
    #import grale.images as images

    #import grale.images as images
    #svd = SceneViewDesciption()
    #l = getPointsLayersFromImagesData(images.ImagesData.load("../../../inversion_examples/example3/images_01.imgdata"), -1)
    #svd.addLayers(l)
    #pprint.pprint(l)

    #img = getImagesDataFromVisibleLayers(svd)
    #pprint.pprint(img)
    #l = getPointsLayersFromImagesData(img, -1)
    #pprint.pprint(l)

#    svd = SceneViewDesciption("/home/jori/projects/a3827new/a3827_tmp.json")
#    for i in range(svd.getNumberOfLayers()):
#        print(i, svd.getLayerName(i))
#
#    for i in [ 3, 5, 6, 7, 8]:
#        svd.setLayerVisibility(i, False)
#    
#    import matplotlib.pyplot as plt
#
#    img = getSceneRegionNumPyArray(svd, [ -20*ANGLE_ARCSEC, -20*ANGLE_ARCSEC ], [ 20*ANGLE_ARCSEC, 20*ANGLE_ARCSEC], 512)
#
#    def showImg(img):
#        plt.figure()
#        plt.imshow(img)
#        plt.gca().invert_yaxis()
#        plt.gca().invert_xaxis()
#        plt.show()
#
#    #showImg(img)
#
#    import pickle
#    ip = pickle.load(open("/home/jori/projects/a3827new/3_pt_ctr_avg-extrasubdiv-nocentral_no.imgplane", "rb"))
#    newLayers, srcAreas = getLayersFromBackProjectRetrace(svd, ip, relensSeparately=False, overWriteFiles=True)
#    pprint.pprint(newLayers)
#
#    for l in newLayers:
#        svd.addLayer(l)
#
#    svd.saveFile("newscene.json")
#
#    img = getSceneRegionNumPyArray(svd, [ -20*ANGLE_ARCSEC, -20*ANGLE_ARCSEC ], [ 20*ANGLE_ARCSEC, 20*ANGLE_ARCSEC], 512)
#    #showImg(img)
#
#    svd = SceneViewDesciption()
#    svd.addLayer(newLayers[0])
#    svd.setAxisVisible(False)

#    img = getSceneRegionNumPyArray(svd, *srcAreas[0], 512)
    #showImg(img)

    svd = SceneViewDesciption("newscene.json")
    src = getSourceShapeFromSourceRGBLayer(svd.getLayer(svd.getNumberOfLayers()-2))

    import grale.lenses as lenses
    import grale.plotutil as plotutil
    from grale.constants import DIST_MPC
    import matplotlib.pyplot as plt

    dummyLens = lenses.PlummerLens(1000*DIST_MPC, { "mass": 1, "width": 100*ANGLE_ARCSEC })

    li = plotutil.LensInfo(dummyLens, size=10*ANGLE_ARCSEC)
    li.setSourceDistances(1,1)
    plotutil.plotImagePlane(li, [src], plotImages=False, sourceRgb=(1,1,1), angularUnit=ANGLE_ARCSEC)
    plt.gca().invert_xaxis()

    plt.show()

if __name__ == "__main__":
    main()
