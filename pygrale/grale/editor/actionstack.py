from __future__ import print_function
import pprint

def objectToDouble(x):
    if type(x) == float:
        return float(x)
    return float("NaN")

class ActionStack(object):
    def __init__(self, maxUndoDepth = 256):
        self.undoStack = [ ]
        self.redoStack = [ ]

        self.maxUndoDepth = maxUndoDepth

    def _getPreviousOperation(self):
        prevOp = None if len(self.undoStack) < 1 else self.undoStack[-1]
        return prevOp

    def _check(self):
        self.redoStack = [ ]
        while len(self.undoStack) >= self.maxUndoDepth:
            self.undoStack.pop(0)

    def clear(self):
        self.undoStack = [ ]
        self.redoStack = [ ]

    def appendToStack(self, x):
        self._check()
        self.undoStack.append(x)

    def recordLabelChange(self, layerUuid, pointUuid, oldLabel, newLabel):
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "label" and prevOp["layer"] == layerUuid and prevOp["point"] == pointUuid:
            prevOp["new"] = newLabel
        else:
            self.appendToStack({ 
                "cmd": "label",
                "layer": layerUuid,
                "point": pointUuid,
                "old": oldLabel,
                "new": newLabel
            })

    def recordPointMove(self, layerUuid, pointUuid, oldPosition, newPosition, moveId):
        pointMove = { 
            "layer": layerUuid,
            "point": pointUuid,
            "old": oldPosition,
            "new": newPosition
        }

        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "move" and moveId and prevOp["id"] == moveId:
            prevOp["moves"].append(pointMove)
        else:
            self.appendToStack({
                "cmd": "move",
                "id": moveId,
                "moves": [ pointMove ]
            })
        

    def recordSetNormalPoint(self, layerUuid, pointUuid, oldPos, newPos, oldTd, newTd, oldLabel, newLabel):
        x = { 
            "cmd": "setpt",
            "layer": layerUuid,
            "point": pointUuid,
            "pos": [oldPos, newPos],
            "timedelay": [oldTd, newTd],
            "label": [oldLabel, newLabel]
        }

        self.appendToStack(x)

    def _addOp(self, addId):
        x = {
            "cmd": "add",
            "id": addId,
            "normalpoints": [ ],
            "matchpoints": [ ],
            "triangles": [ ]
        }
        self.appendToStack(x)
        return x

    def recordAddNormalPoint(self, layerUuid, pointUuid, pos, timedelay, label, addId):
        pointAdd = { 
            "layer": layerUuid,
            "point": pointUuid,
            "pos": pos,
            "timedelay": timedelay,
            "label": label
        }

        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "add" and addId and prevOp["id"] == addId:
            prevOp["normalpoints"].append(pointAdd)
        else:
            self._addOp(addId)["normalpoints"].append(pointAdd)

    def recordAddNormalPoints(self, points, addId):
        x = {
            "cmd": "add",
            "id": addId,
            "normalpoints": points,
            "matchpoints": [ ],
            "triangles": [ ]
        }
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "add" and addId and prevOp["id"] == addId:
            prevOp["normalpoints"] += points
        else:
            self.appendToStack(x)
        return x

    def recordAddMatchPoints(self, points, addId):
        x = {
            "cmd": "add",
            "id": addId,
            "normalpoints": [ ],
            "matchpoints": points,
            "triangles": [ ]
        }
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "add" and addId and prevOp["id"] == addId:
            prevOp["matchpoints"] += points
        else:
            self.appendToStack(x)
        return x

    def recordAddTriangle(self, layerUuid, triangUuid, pointUuids, addId):
        triangAdd = {
            "triangle": triangUuid,
            "layer": layerUuid,
            "points": pointUuids
        }

        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "add" and addId and prevOp["id"] == addId:
            prevOp["triangles"].append(triangAdd)
        else:
            self._addOp(addId)["triangles"].append(triangAdd)

    def recordAddTriangles(self, triangles, addId):
        x = {
            "cmd": "add",
            "id": addId,
            "normalpoints": [ ],
            "matchpoints": [ ],
            "triangles": triangles
        }
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "add" and addId and prevOp["id"] == addId:
            prevOp["triangles"] += triangles
        else:
            self.appendToStack(x)
        return x

    def recordAddMatchPoint(self, layerUuid, pointUuid, pos, label, addId):
        matchAdd = {
            "layer": layerUuid,
            "point": pointUuid,
            "pos": pos,
            "label": label
        }

        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "add" and addId and prevOp["id"] == addId:
            prevOp["matchpoints"].append(matchAdd)
        else:
            self._addOp(addId)["matchpoints"].append(matchAdd)

    def _delOp(self, deleteId):
        x = {
            "cmd": "del",
            "id": deleteId,
            "normalpoints": [ ],
            "matchpoints": [ ],
            "triangles": [ ]
        }
        self.appendToStack(x)
        return x

    def recordDeletePoints(self, points, deleteId):
        x = {
            "cmd": "del",
            "id": deleteId,
            "normalpoints": points,
            "matchpoints": [ ],
            "triangles": [ ]
        }
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "del" and deleteId and prevOp["id"] == deleteId:
            prevOp["normalpoints"] += points
        else:
            self.appendToStack(x)
        return x

    def recordDeleteMatchPoints(self, points, deleteId):
        x = {
            "cmd": "del",
            "id": deleteId,
            "normalpoints": [ ],
            "matchpoints": points,
            "triangles": [ ]
        }
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "del" and deleteId and prevOp["id"] == deleteId:
            prevOp["matchpoints"] += points
        else:
            self.appendToStack(x)
            pprint.pprint(x)
        return x

    def recordDeleteTriangles(self, triangs, deleteId):
        x = {
            "cmd": "del",
            "id": deleteId,
            "normalpoints": [ ],
            "matchpoints": [ ],
            "triangles": triangs
        }
        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "del" and deleteId and prevOp["id"] == deleteId:
            prevOp["triangles"] += triangs
        else:
            self.appendToStack(x)
        return x

    def recordDeletePoint(self, layerUuid, pointUuid, pos, timedelay, label, deleteId):
        x = { 
            "layer": layerUuid,
            "point": pointUuid,
            "pos": pos,
            "timedelay": timedelay,
            "label": label
        }

        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "del" and deleteId and prevOp["id"] == deleteId:
            prevOp["normalpoints"].append(x)
        else:
            self._delOp(deleteId)["normalpoints"].append(x)

    def recordDeleteTriangle(self, layerUuid, triangUuid, pointUuids, deleteId):
        x = { 
            "layer": layerUuid,
            "triangle": triangUuid,
            "points": pointUuids,
        }

        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "del" and deleteId and prevOp["id"] == deleteId:
            prevOp["triangles"].append(x)
        else:
            self._delOp(deleteId)["triangles"].append(x)

    def recordDeleteMatchPoint(self, layerUuid, pointUuid, pos, label, deleteId):
        x = { 
            "layer": layerUuid,
            "point": pointUuid,
            "pos": pos,
            "label": label
        }

        # Check if we can merge
        prevOp = self._getPreviousOperation()
        if prevOp and not "finished" in prevOp and prevOp["cmd"] == "del" and deleteId and prevOp["id"] == deleteId:
            prevOp["matchpoints"].append(x)
        else:
            self._delOp(deleteId)["matchpoints"].append(x)

    def recordCenterLocation(self, layerUuid, oldCenter, newCenter):
        self.appendToStack({
            "cmd": "center",
            "layer": layerUuid,
            "old": oldCenter,
            "new": newCenter
        })

    def recordCenterAndMinMax(self, layerUuid, oldCenter, newCenter, oldMinMax, newMinMax):
        self.appendToStack({
            "cmd": "centerminmax",
            "layer": layerUuid,
            "oldcenter": oldCenter,
            "newcenter": newCenter,
            "oldminmax": oldMinMax,
            "newminmax": newMinMax
        })

    def recordRGBTransform(self, layerUuid, oldTransform, newTransform):
        self.appendToStack({
            "cmd": "transform",
            "layer": layerUuid,
            "old": oldTransform,
            "new": newTransform
        })

    def _undoRedoLabel(self, isRedo, x, scene):
        layerItem = scene.getLayerItem(x["layer"])
        label = x["new"] if isRedo else x["old"]
        layerItem.setLabel(x["point"], label)

    def _undoRedoMove(self, isRedo, x, scene):
        for pointMove in x["moves"]:
            layerUuid = pointMove["layer"]
            uuid = pointMove["point"]
            pos = pointMove["new"] if isRedo else pointMove["old"]
            layerItem = scene.getLayerItem(layerUuid)
            layerItem.movePoint(uuid, pos[0], pos[1])

    def _redoAdd(self, x, scene):
        for p in x["normalpoints"]:
            layerUuid = p["layer"]
            layerItem = scene.getLayerItem(layerUuid)
            layerItem.addPoint(p["pos"][0], p["pos"][1], p["label"], objectToDouble(p["timedelay"]), p["point"])

        for t in x["triangles"]:
            layerUuid = t["layer"]
            layerItem = scene.getLayerItem(layerUuid)
            layerItem.addTriangle(*t["points"], t["triangle"])

        for m in x["matchpoints"]:
            layerUuid = m["layer"]
            layerItem = scene.getLayerItem(layerUuid)
            layerItem.addPoint(m["pos"][0], m["pos"][1], m["label"], objectToDouble(None), m["point"])

    def _redoDelete(self, x, scene):

        for t in x["triangles"]:
            layerUuid = t["layer"]
            layerItem = scene.getLayerItem(layerUuid)
            layerItem.clearTriangle(t["triangle"])

        for p in x["normalpoints"]:
            layerUuid = p["layer"]
            layerItem = scene.getLayerItem(layerUuid)
            layerItem.clearPoint(p["point"])

        for m in x["matchpoints"]:
            layerUuid = m["layer"]
            layerItem = scene.getLayerItem(layerUuid)
            layerItem.clearPoint(m["point"])

    def _undoRedoCenter(self, isRedo, x, scene):
        
        ra, dec = x["new"] if isRedo else x["old"]
        layerItem = scene.getLayerItem(x["layer"])
        layer = layerItem.getLayer()
        layer.setCenter(ra, dec)
        layerItem.updateTransform()

    def _undoRedoCenterMinMax(self, isRedo, x, scene):
        
        ra, dec = x["newcenter"] if isRedo else x["oldcenter"]
        minValue, maxValue = x["newminmax"] if isRedo else x["oldminmax"]
        layerItem = scene.getLayerItem(x["layer"])
        layer = layerItem.getLayer()
        layer.setCenter(ra, dec)
        layer.setMinMax(minValue, maxValue)
        layerItem.updatePixmap()
        layerItem.updateTransform()

    def _undoRedoTransform(self, isRedo, x, scene):

        t = x["new"] if isRedo else x["old"]
        layerItem = scene.getLayerItem(x["layer"])
        layer = layerItem.getLayer()
        layer.setImageTransform(t)
        layerItem.updateTransform()

    def _undoRedoSetPoint(self, isRedo, x, scene):
        pos = x["pos"][1] if isRedo else x["pos"][0]
        label = x["label"][1] if isRedo else x["label"][0]
        timedelay = x["timedelay"][1] if isRedo else x["timedelay"][0]
        layerItem = scene.getLayerItem(x["layer"])
        layerItem.setPoint(x["point"], pos[0], pos[1], label, objectToDouble(timedelay))

    def _undoRedo(self, isRedo, x, scene):
        if x["cmd"] == "label":
            self._undoRedoLabel(isRedo, x, scene)
        elif x["cmd"] == "move":
            self._undoRedoMove(isRedo, x, scene)
        elif x["cmd"] == "add":
            if isRedo:
                self._redoAdd(x, scene)
            else:
                self._redoDelete(x, scene)
        elif x["cmd"] == "del":
            if isRedo:
                self._redoDelete(x, scene)
            else:
                self._redoAdd(x, scene)
        elif x["cmd"] == "center":
            self._undoRedoCenter(isRedo, x, scene)
        elif x["cmd"] == "transform":
            self._undoRedoTransform(isRedo, x, scene)
        elif x["cmd"] == "setpt":
            self._undoRedoSetPoint(isRedo, x, scene)
        elif x["cmd"] == "centerminmax":
            self._undoRedoCenterMinMax(isRedo, x, scene)
        else:
            self.customStackEntryUndoRedo(x, isRedo, scene)

    def undo(self, scene):
        if len(self.undoStack) == 0:
            return

        # Move the item from the undo list to the redo list
        x = self.undoStack.pop()
        self.redoStack.insert(0, x)
        x["finished"] = True

        self._undoRedo(False, x, scene)

    def redo(self, scene):
        if len(self.redoStack) == 0:
            return

        x = self.redoStack.pop(0)
        self.undoStack.append(x)
        x["finished"] = True

        self._undoRedo(True, x, scene)

    def customStackEntryUndoRedo(self, x, isRedo, scene):
        raise Exception("Unknown command: {}".format(x))

