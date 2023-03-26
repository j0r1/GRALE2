from base import GraphicsView, Layer
import scenes
import sys
import json
import pprint
import multiprocessing
import checkqt
import tools
from multiprocessing import Process, Queue
import copy
import os

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

    def addLayers(self, layers, visibilities):
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
            for k in [ "hduidx", "center", "minma" ]:
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
    checkqt.checkQtAvailable()

    q1, q2 = Queue(), Queue()
    p = Process(target=helperFunction, args=(q1,q2))
    p.start()
    q1.put(paramList)
    flag, result = q2.get()
    p.join()
    if not flag:
        raise Exception(result)
    return result

# Coords are in arcsec for now
def getSceneRegionNumPyArray(sceneViewDesc, bottomLeft, topRight, widthPixels = None, heightPixels = None, grayScale = False):
    return _foreground_helper(_background_helper_getSceneRegionImageNumPy,
                              [ sceneViewDesc.toObject(), bottomLeft, topRight, widthPixels, heightPixels, grayScale ])

def _background_helper_getPointsLayersFromImagesData(q1, q2):
    _background_helper(q1, q2, lambda params: [ l.toSettings() for l in tools.importImagesDataToLayers(*params) ])

def getPointsLayersFromImagesData(imgDat, which, layerTitleFormat = "Points from images data"):
    return _foreground_helper(_background_helper_getPointsLayersFromImagesData, [ imgDat, which, layerTitleFormat ])

def main():

    #svd = SceneViewDesciption("cl0024.json")
    #pprint.pprint(svd.toObject())

    #svd = SceneViewDesciption("cl0024.json")
    #img = getSceneRegionNumPyArray(svd, [ -10, -10 ], [10, 10], 512)
    #
    #import matplotlib.pyplot as plt
    #plt.imshow(img)
    #plt.show()
    import grale.images as images

    l = getPointsLayersFromImagesData(images.ImagesData.load("../../../inversion_examples/example3/images_01.imgdata"), -1)
    pprint.pprint(l)

if __name__ == "__main__":
    main()
