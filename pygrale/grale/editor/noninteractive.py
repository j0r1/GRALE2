from base import GraphicsView, Layer
import scenes
import sys
import json
import pprint
import multiprocessing
import checkqt
from multiprocessing import Process, Queue

class MultipleLayerScene(scenes.LayerScene):

    def __init__(self):
        super(MultipleLayerScene, self).__init__()

        self.layersAndWidgets = []
        self.nextZ = 0

    def addLayer(self, l):
        w = l.createGraphicsItem()
        self.layersAndWidgets.append((l, w))

        self.nextZ += 1
        w.setZValue(self.nextZ)
        return w

    def getLayerItem(self, uuid):
        if uuid in self.layerItems:
            return self.layerItems[uuid]
        raise Exception("Specified layer {} not found".format(uuid))

    def getCurrentItemAndLayer(self):
        return self.layersAndWidgets[-1]

def importLayers(scene, layers, visibilities):
    if visibilities is None:
        visibilities = [ True for _ in layers ]

    for layer,visible in zip(layers,visibilities):
        if not visible:
            continue

        l = Layer.fromSettings(layer)
        #pprint.pprint(l)
        w = scene.addLayer(l)
        scene.addItem(w)
        w.updatePoints()

def setGlobalSettings(scene, d):
    scene.setAxesVisible(d["showaxes"])
    scene.setAxisLeft(d["axisleft"])
    scene.setPointSizePixelSize(d["pointsizepixelsize"])
    scene.setPointSizeSceneSize(d["pointsizearcsecsize"])
    scene.setPointSizeFixed(d["ispointsizefixed"])
    scene.setMatchPointsVisible(d["showmatchpoints"])

def setViewSettings(view, v):
    ctr = v["center"]
    s = v["zoom"]
    view.centerOn(ctr[0], ctr[1])
    view.setScale(s)

def _loadFile(fileName):
    from PyQt5 import QtWidgets, QtCore, QtGui

    app = QtWidgets.QApplication(sys.argv)
    scene = MultipleLayerScene()
    #view = GraphicsView(scene)
    
    obj = json.load(open(fileName, "rt"))
    if "state" in obj:
        layers = obj["state"]["layers"]
        visibilities = obj["state"]["visibilities"]
        setGlobalSettings(scene, obj["state"]["globalsettings"])
    else:
        layers = obj
        visibilities = None

    scale = 1 if not "view" in obj else obj["view"]["zoom"]
    scene.onScaleChanged(scale) # May be needed to recalculate the point size transformation

    importLayers(scene, layers, visibilities)

    return app, scene

def _helper_getSceneRegionImageNumPy(q1, q2):
    fileName, bottomLeft, topRight, widthPixels, heightPixels, grayScale = q1.get()

    try:
        app, scene = _loadFile(fileName)
        res = scene.getSceneRegionNumPyArray(bottomLeft, topRight, widthPixels, heightPixels, grayScale)
        #print(res)
        q2.put([True, res])
    except Exception as e:
        q2.put([False, str(e)])

# Coords are in arcsec for now
def getSceneRegionNumPyArray(fileName, bottomLeft, topRight, widthPixels = None, heightPixels = None, grayScale = False):
    checkqt.checkQtAvailable()

    q1, q2 = Queue(), Queue()
    p = Process(target=_helper_getSceneRegionImageNumPy, args=(q1,q2))
    p.start()
    q1.put([fileName, bottomLeft, topRight, widthPixels, heightPixels, grayScale ])
    flag, result = q2.get()
    p.join()
    if not flag:
        raise Exception(result)
    return result
    
def main():

    img = getSceneRegionNumPyArray("cl0024.json", [ -10, -10 ], [10, 10], 512)

    import matplotlib.pyplot as plt
    plt.imshow(img)
    plt.show()

if __name__ == "__main__":
    main()
