from fractions import Fraction
from . import gridfunction
import pprint
import numpy as np
import sys
import copy

class GridException(Exception):
    pass

_defaultExcludeFunction = lambda pos, size: False

def _realCellCenterAndSize(size, center, cell):
    cellSize = float(cell["size"])*size
    x = center[0] + float(cell["center"][0])*size
    y = center[1] + float(cell["center"][1])*size
    return ([x, y], cellSize)

def _fractionalGridToRealGrid(grid):
    g = [ ]
    size = grid["size"]
    center = grid["center"]
    for cell in grid["cells"]:
        realPos, realSize = _realCellCenterAndSize(size, center, cell)
        g.append({ "center": realPos, "size": realSize })
    return g

def createUniformGrid(size, center, axisSubDivisions, excludeFunction = _defaultExcludeFunction):

    if axisSubDivisions < 1:
        raise GridException("The number of axis subdivisions for a uniform grid must be at least 1")

    grid = { "size": size, "center": center, "cells": [ ] }
    doubleCellSize = Fraction(2, axisSubDivisions)
    cellSize = doubleCellSize/2
    for y in range(axisSubDivisions):
        # This calculates from -1 to 1
        cellCenterY = Fraction(-1,1) + Fraction(1, axisSubDivisions) + y*doubleCellSize
        cellCenterY /= 2 # to make it go from -1/2 to 1/2 (so that the width is 1)

        for x in range(axisSubDivisions):
            cellCenterX = Fraction(-1,1) + Fraction(1, axisSubDivisions) + x*doubleCellSize
            cellCenterX /= 2

            cell = {
                "center": [ cellCenterX, cellCenterY ],
                "size": cellSize
            }

            realPos, realSize = _realCellCenterAndSize(size, center, cell)
            if not excludeFunction(realPos, realSize):
                grid["cells"].append(cell)

    return grid

def _integrateGridFunction(f, center, size, minPointDist):
    num = int(size/minPointDist)
    if num < 2:
        num = 2

    bottomLeft = np.array([center[0]-size/2.0, center[1]-size/2.0])
    dx = size/float(num-1)
    dy = dx
    dxy = np.array([ dx, dy ])

    #points = np.zeros([num, num, 2])
    #for y in range(num):
    #    for x in range(num):
    #        points[y,x,0] = bottomLeft[0]+x*dxy[0]
    #        points[y,x,1] = bottomLeft[0]+y*dxy[1]

    # Numpy way to do the above, not sure if it can be done in a better way
    points = np.empty([num, num, 2])
    xcoords = np.linspace(center[0]-size/2.0, center[0]+size/2.0, num)
    ycoords = np.linspace(center[1]-size/2.0, center[1]+size/2.0, num)
    points[:,:,0] = np.outer(np.ones(ycoords.size), xcoords)
    points[:,:,1] = np.outer(ycoords, np.ones(xcoords.size))

    densities = f.evaluate(points)
    pixels = 0.25*(densities[0:num-1,0:num-1] + densities[0:num-1,1:num] + 
                   densities[1:num,0:num-1] + densities[1:num,1:num])
    #pprint.pprint(densities)
    #pprint.pprint(pixels)

    pixels *= dxy[0]*dxy[1]
    return pixels.sum()

# TODO: use excludeFunction?
def _createSubdivisionGridForThreshold(f, size, center, thresholdMass, startSubDiv, excludeFunction,
        maxIntegrationSubDiv, keepLarger, cache):
    
    grid = createUniformGrid(size, center, startSubDiv, excludeFunction)
    gridCells = grid["cells"]
    minPointDist = size/float(maxIntegrationSubDiv)

    count = 0
    maxCount = 4096
    while count < maxCount:
        count += 1
        #print("count = ", count)
        newCells = [ ]
        gotSubDiv = False

        for cell in gridCells:
            if "marked" in cell: 
                # just keeping a larger cell that's been subdivided, or one that did not
                # contain enough mass
                newCells.append(cell)
            else:
                # Check if we need to subdivide this
                realPos, realSize = _realCellCenterAndSize(size, center, cell)
                key = (realPos[0], realPos[1], realSize)
                if key in cache:
                    cellMass = cache[key]
                else:
                    cellMass = _integrateGridFunction(f, realPos, realSize, minPointDist)
                    cache[key] = cellMass

                #print("cellMass=",cellMass, "thresholdMass=",thresholdMass)
                if cellMass > thresholdMass:

                    gotSubDiv = True

                    if keepLarger:
                        cell["marked"] = True
                        newCells.append(cell)

                    origSize = cell["size"]
                    origX, origY = cell["center"]
                    newCells.append({ "size": origSize/2,
                                      "center": [ origX-origSize/4, origY-origSize/4 ] })
                    newCells.append({ "size": origSize/2,
                                      "center": [ origX-origSize/4, origY+origSize/4 ] })
                    newCells.append({ "size": origSize/2,
                                      "center": [ origX+origSize/4, origY+origSize/4 ] })
                    newCells.append({ "size": origSize/2,
                                      "center": [ origX+origSize/4, origY-origSize/4 ] })

                else:
                    cell["marked"] = True
                    newCells.append(cell)

        gridCells = newCells
        if not gotSubDiv:
            # We're done, clean up the 'marked' flags
            for cell in gridCells:
                if "marked" in cell:
                    del cell["marked"]

            grid["cells"] = gridCells
            return grid

    raise GridException("Unexpected: couldn't find a subdivision grid within {} iterations".format(maxCount))

def createSubdivisionGridForFITS(fitsHDUEntry, centerRaDec, gridSize, gridCenter, minSquares, maxSquares, startSubDiv = 1,
        excludeFunction = _defaultExcludeFunction,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        ignoreOffset = True,
        useAbsoluteValues = True):
    
    class tmpClass(object):
        pass

    hduCopy = tmpClass
    hduCopy.data = fitsHDUEntry.data.copy()
    hduCopy.header = fitsHDUEntry.header[:]

    data = hduCopy.data
    if useAbsoluteValues:
        data = np.absolute(data)

    minValue = data.min()
    if minValue < 0:
        print("Warning: negative values will be clipped")
        data = data.clip(0, None)
        minValue = 0

    if ignoreOffset:
        data = data - minValue

    f = gridfunction.GridFunction.createFromFITS(hduCopy, centerRaDec, True)
    return _createSubdivisionGridCommon(f, gridSize, gridCenter, minSquares, maxSquares, startSubDiv, 
                                        excludeFunction, maxIntegrationSubDiv, keepLarger)

def _createSubdivisionGridCommon(f, size, center, minSquares, maxSquares, startSubDiv = 1, 
        excludeFunction = _defaultExcludeFunction,
        maxIntegrationSubDiv = 256,
        keepLarger = False):

    if minSquares >= maxSquares:
        raise GridException("Minimal number of cells must be smaller than maximum")
    if maxSquares < startSubDiv**2:
        raise GridException("Requested initial subdivision leads to more cells than the specified maximum")

    # Note that this isn't really the total mass, we should multiply with Dd**2 for this
    minPointDist = size/float(maxIntegrationSubDiv)
    totalMass = _integrateGridFunction(f, center, size, minPointDist)
    #print("total mass in area:", totalMass*lensInfo["Dd"]**2)

    subDivFraction = 0.2
    diffFrac = 2.0
    maxDiffFracIt = 50
    maxSubDivIt = 1000

    massCache = { }

    for i in range(maxDiffFracIt):
        diffFrac /= 2.0

        maxIt = maxSubDivIt
        grid = None

        while not grid or (len(grid["cells"]) < minSquares and maxIt > 0):
            maxIt -= 1
            subDivFraction /= (1.0 + diffFrac)
            grid = _createSubdivisionGridForThreshold(f, size, center, totalMass*subDivFraction, startSubDiv, 
                                                      excludeFunction, maxIntegrationSubDiv, keepLarger, massCache)

        if len(grid["cells"]) < minSquares:
            raise GridException("Unable to find a grid that has a larger number of squares than the specified minimum")

        if len(grid["cells"]) <= maxSquares:
            return grid

        diffFrac /= 2.0
        maxIt = maxSubDivIt

        while len(grid["cells"]) > maxSquares and maxIt > 0:
            maxIt -= 1
            subDivFraction *= (1.0 + diffFrac)
            grid = _createSubdivisionGridForThreshold(f, size, center, totalMass*subDivFraction, startSubDiv, 
                                                      excludeFunction, maxIntegrationSubDiv, keepLarger, massCache)

        if len(grid["cells"]) > maxSquares:
            raise GridException("Unable to find a grid that has a smaller number of squares than the specified maximum")

        if len(grid["cells"]) >= minSquares:
            return grid

    raise GridException("A grid of the requested size could not be constructed")

def createSubdivisionGrid(size, center, lensInfo, minSquares, maxSquares, startSubDiv = 1, 
        excludeFunction = _defaultExcludeFunction,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        ignoreOffset = True,
        useAbsoluteValues = True):

    densPoints = lensInfo["densitypoints"]
    if useAbsoluteValues:
        densPoints = np.absolute(densPoints)

    offset = densPoints.min()
    if offset < 0:
        raise GridException("Detected negative value in density points, can't handle this")

    if ignoreOffset:
        densPoints = densPoints - offset

    f = gridfunction.GridFunction(densPoints, lensInfo["bottomleft"], lensInfo["topright"])

    return _createSubdivisionGridCommon(f, size, center, minSquares, maxSquares, startSubDiv, 
                                        excludeFunction, maxIntegrationSubDiv, keepLarger)

def fitMultiplePlummerLens(gridCells, Dd, targetDensityFunction, sizeFactor = 1.7):

    from . import lenses

    # Use a set, and fractional calculus (if available) to avoid double points
    points = set([])
    for cell in gridCells["cells"]:
        cellSize = cell["size"]
        cellSize2 = cellSize/2
        center = cell["center"]
        points.add((center[0], center[1]))
        points.add((center[0] + cellSize2, center[1] ))
        points.add((center[0] - cellSize2, center[1] ))
        points.add((center[0], center[1] + cellSize2 ))
        points.add((center[0], center[1] - cellSize2 ))

    # Convert to real values
    gridSize = float(gridCells["size"])
    gridCenter = np.array([float(gridCells["center"][0]), float(gridCells["center"][1])])
    points = np.array([ np.array([float(p[0]), float(p[1])]) * gridSize + gridCenter for p in points])

    cells = _fractionalGridToRealGrid(gridCells)
    plummers = [ ]
    for cell in cells:
        w = cell["size"]*sizeFactor
        m = np.pi * Dd**2*w**2
        plum = lenses.PlummerLens(Dd, { "mass": m, "width": w })
        plummers.append([ np.array(cell["center"]), plum ])
    
    def evaluatePlummers(plummers, points):
        B = np.empty(shape=[len(plummers),len(points)], dtype=np.double)
        for i in range(len(plummers)):
            plum = plummers[i]
            center = plum[0]
            basisLens = plum[1]
            relPoints = points - center
            B[i,:] = basisLens.getSurfaceMassDensity(relPoints)

        return np.asmatrix(B)
    
    B = evaluatePlummers(plummers, points)
    f = targetDensityFunction(points)

    a = B.dot(B.transpose())
    b = B.dot(f).reshape((-1,1))
    
    weights = np.linalg.solve(a, b)
    #print(a.shape, b.shape, f.shape)

    multiplePlummerParams = [ ]
    for i in range(len(plummers)):
        plum = plummers[i]
        center = plum[0]
        basisLens = plum[1]
        params = basisLens.getLensParameters()
        m = params["mass"] * weights[i]
        w = params["width"]
        multiplePlummerParams.append({
            "mass": m,
            "width": w,
            "x": center[0],
            "y": center[1]
        })

    fitLens = lenses.MultiplePlummerLens(Dd, multiplePlummerParams)
    
    #print("f.shape", f.shape)
    #print("B.transpose().shape", B.transpose().shape)
    #print("weights.shape", weights.shape)
    #predDiff = (f.reshape((-1,1)) - B.transpose().dot(weights))
    #print("predDiff.shape", predDiff.shape)
    #s = predDiff.transpose().dot(predDiff)
    #x2 = float(s)/len(points)

    #return (fitLens, x2)
    return fitLens

def _debugPlot(g):
    for cell in g:
        size = cell["size"]
        x, y = cell["center"]
        print(x-size/2, y-size/2)
        print(x-size/2, y+size/2)
        print(x+size/2, y+size/2)
        print(x+size/2, y-size/2)
        print(x-size/2, y-size/2)
        print()

def main():
    #grid = createUniformGrid(123, [-1, -2], 5)
    #pprint.pprint(grid)
    #g = _fractionalGridToRealGrid(grid)
    #_debugPlot(g)

    import grale.lenses as lenses
    import grale.plotutil as plotutil
    from grale.constants import ANGLE_ARCSEC

    U = ANGLE_ARCSEC
    lensInfo = {
        "lens": lenses.GravitationalLens.load("/home/jori/projects/gamods-hg/examples/reallens_nosheet.lensdata"),
        "bottomleft": [-100*U, -100*U],
        "topright": [100*U, 100*U],
    }
    plotutil.plotDensity(lensInfo, feedbackObject = "none")
    #pprint.pprint(lensInfo)

    grid = createSubdivisionGrid(190*U, [0,0], lensInfo, 900, 1000, startSubDiv=7)
    print("# num cells", len(grid["cells"]))
    #pprint.pprint(grid)

    g = _fractionalGridToRealGrid(grid)
    _debugPlot(g)

if __name__ == "__main__":
    main()
