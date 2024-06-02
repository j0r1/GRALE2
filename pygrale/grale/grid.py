"""This module contains functions to create a subdivision grid that's used
in the lens inversion procedure. Apart from a uniformly subdivided grid, it
is also possible to lay out the grid based on some kind of density map, using 
a finer subdivision where the integrated density is higher. This density
could be the projected mass density of a gravitational lens, or even the
observed light distribution in a FITS file.

Examples for how these functions can be used, can be found in the notebooks
`gridtests.ipynb <_static/gridtests.ipynb>`_ and
`fittest.ipynb <_static/fittest.ipynb>`_.
"""
from fractions import Fraction
from . import gridfunction
import numpy as np
import random
import copy

class GridException(Exception):
    """An exception that's generated when something goes wrong with
    the grid creation."""
    pass

_defaultExcludeFunction = lambda pos, size, inUniformStep: False
_defaultSubDivFunction = lambda pos, size: True

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

def fractionalGridToRealGrid(grid):
    """Transforms the specified `grid`, which uses a grid center and
    fractional sizes and centers for the cells, to a straighforward list of 
    cells, each having a regular center and size.
    """ 
    return _fractionalGridToRealGrid(grid)

def createMultiUniformGrid(uniformRegionInfo, # [ { size, center, axissubdiv } ]
                           excludeFunction = None):
    """TODO"""
    
    if not uniformRegionInfo:
        raise GridException("No regions were specified")

    for inf in uniformRegionInfo:
        if inf["axissubdiv"] < 1:
            raise GridException("The number of axis subdivisions for each uniform grid must be at least 1")

    if not excludeFunction:
        excludeFunction = _defaultExcludeFunction

    allGrids = []
    for inf in uniformRegionInfo:
        size, center, axisSubDivisions = inf["size"], inf["center"], inf["axissubdiv"]

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
                if not excludeFunction(realPos, realSize, True):
                    grid["cells"].append(cell)

        allGrids.append(grid)

    return allGrids

def createUniformGrid(size, center, axisSubDivisions, excludeFunction = None):
    """Creates a uniform grid, centered on the X,Y coordinates in `center` and
    with width and height specified by `size`. The generated grid will be have
    `axisSubDivisions` cells in X and Y directions.

    If `excludeFunction` is specified, for each grid cell it will be called as
    ``excludeFunction(cellCenter, cellSize, True)``, and if the function returns ``True``
    then that particular cell will not be included in the final grid.
    """

    return createMultiUniformGrid([{ "center": center, "size": size, "axissubdiv": axisSubDivisions}], 
                                  excludeFunction)[0]

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

    densities = f(points)
    pixels = 0.25*(densities[0:num-1,0:num-1] + densities[0:num-1,1:num] + 
                   densities[1:num,0:num-1] + densities[1:num,1:num])
    #pprint.pprint(densities)
    #pprint.pprint(pixels)

    pixels *= dxy[0]*dxy[1]
    return pixels.sum()

def _getNewSubDivCells(f, gridCells, size, center, minPointDist, cache, thresholdMass, keepLarger,
                       checkSubDivFunction):
        
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

            shouldSubDiv = False
            if checkSubDivFunction(realPos, realSize): # We're allowed to subdivide this

                key = (realPos[0], realPos[1], realSize)
                if key in cache:
                    cellMass = cache[key]
                else:
                    cellMass = _integrateGridFunction(f, realPos, realSize, minPointDist)
                    cache[key] = cellMass

                if cellMass > thresholdMass:
                    shouldSubDiv = True

            #print("cellMass=",cellMass, "thresholdMass=",thresholdMass)
            if shouldSubDiv:

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

    return newCells, gotSubDiv

def _checkGridExcludeFunction(grid, gridCells, size, center, excludeFunction):
    grid["cells"] = [ ]
    for cell in gridCells:
        if "marked" in cell:
            del cell["marked"]

        realPos, realSize = _realCellCenterAndSize(size, center, cell)
    
        if not excludeFunction(realPos, realSize, False):
            grid["cells"].append(cell)
    
def _createSubdivisionGridForThreshold_Multi(f, 
        # { "size": , "center": , "startsubdiv": }
        subdivRegionInfo, # list of these ^ things
        thresholdMass,
        excludeFunction,
        maxIntegrationSubDiv, keepLarger, cache, checkSubDivFunction, maxSquares):
    
    if not subdivRegionInfo:
        raise GridException("No subdivision regions were specified")

    allGrids = [ createUniformGrid(x["size"], x["center"], x["startsubdiv"], excludeFunction) for x in subdivRegionInfo ]
    allGridCells = [ copy.copy(grid["cells"]) for grid in allGrids ]
    allGridCellCount = sum([ len(cells) for cells in allGridCells ])
    
    assert type(f) == list, "Expecting a list of functions"
    if len(f) != len(subdivRegionInfo):
        raise GridException("Expecting as many functions as regions")
    allFunctions = f # Rename so we can reuse f below

    if allGridCellCount > maxSquares:
        raise GridException("Initial uniform subdivision leads to more cells than the specified maximum")

    minPointDists = [ x["size"]/float(maxIntegrationSubDiv) for x in subdivRegionInfo ]

    count = 0
    maxCount = 4096
    while count < maxCount:
        count += 1
        #print("count = ", count)
        allNewCells = [ ]
        gotSubDiv = False

        for f, gridCells, minPointDist, regionInfo in zip(allFunctions, allGridCells, minPointDists, subdivRegionInfo):

            newCells, gotSubDivPart = _getNewSubDivCells(f, gridCells, regionInfo["size"], regionInfo["center"],
                                          minPointDist, cache, thresholdMass, keepLarger,
                                          checkSubDivFunction)
            if gotSubDivPart:
                gotSubDiv = True

            allNewCells.append(newCells)
        
        allGridCells = allNewCells
        if not gotSubDiv:
            
            # We're done, check the excludeFunction

            for f, grid, gridCells, regionInfo in zip(allFunctions, allGrids, allGridCells, subdivRegionInfo):
            
                _checkGridExcludeFunction(grid, gridCells, regionInfo["size"],
                                          regionInfo["center"], excludeFunction)            
            return allGrids

    raise GridException("Unexpected: couldn't find a subdivision grid within {} iterations".format(maxCount))

def createMultiSubdivisionGridForFITS(fitsHDUEntry, centerRaDec,
        subdivRegionInfo, # [ { size, center, startsubdiv } ]
        minSquares, 
        maxSquares,
        excludeFunction = None,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        ignoreOffset = True,
        useAbsoluteValues = True,
        checkSubDivFunction = None):
    """TODO"""

    class tmpClass(object):
        pass

    if not excludeFunction:
        excludeFunction = _defaultExcludeFunction
    if not checkSubDivFunction:
        checkSubDivFunction = _defaultSubDivFunction

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
    return createMultiSubdivisionGridForFunction(f.evaluate, subdivRegionInfo, minSquares, maxSquares,
                                        excludeFunction, maxIntegrationSubDiv, keepLarger, checkSubDivFunction)

def createSubdivisionGridForFITS(fitsHDUEntry, centerRaDec, gridSize, gridCenter, minSquares, maxSquares, startSubDiv = 1,
        excludeFunction = None,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        ignoreOffset = True,
        useAbsoluteValues = True,
        checkSubDivFunction = None):
    """Creates a subdivision grid that's based on an entry of a FITS file. It calls
    :func:`createSubdivisionGridForFunction` internally, where you can find the
    explanation of the algorithm used.

    Arguments:

     - `fitsHDUEntry`: the FITS entry on which the subdivision grid should be based
     - `centerRaDec`: use these RA,Dec coordinates to recalculate the coordinates in the
       FITS file, i.e. that specific point will be associated to the new (0,0) coordinate.
     - `gridSize`: the width and height of the grid, will be used as the `size` parameter
       in :func:`createSubdivisionGridForFunction`
     - `gridCenter`: the center of the grid in the re-centered coordinates. Will be used
       as the `center` parameter in :func:`createSubdivisionGridForFunction`.
     - `minSquares`, `maxSquares`, `startSubDiv`, `excludeFunction`, `maxIntegrationSubDiv`,
       `keepLarger`: see :func:`createSubdivisionGridForFunction`
     - `ignoreOffset`: adds or subtracts a value from the data so that the minimum value
       becomes zero.
     - `useAbsoluteValues`: set to ``False`` to allow negative values (probably not a good
       idea).
     - `checkSubDivFunction`: TODO
    """

    return createMultiSubdivisionGridForFITS(fitsHDUEntry, centerRaDec,
                                             [{ "size": gridSize, "center": gridCenter, "startsubdiv": startSubDiv }],
                                             minSquares, maxSquares, excludeFunction, maxIntegrationSubDiv,
                                             keepLarger, ignoreOffset, useAbsoluteValues, checkSubDivFunction)[0]

def createMultiSubdivisionGridForFunction(targetDensityFunction,
        subdivRegionInfo, # [ { size, center, startsubdiv } ]
        minSquares, maxSquares,
        excludeFunction = None,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        checkSubDivFunction = None):
    """TODO"""

    if not subdivRegionInfo:
        raise GridException("No subdivision region info was specified")
    
    subdivRegionInfo = copy.copy(subdivRegionInfo)
    for x in subdivRegionInfo:
        if not "startsubdiv" in x:
            x["startsubdiv"] = 1

    if type(targetDensityFunction) != list:
        targetDensityFunction = [ targetDensityFunction for x in subdivRegionInfo ]

    if len(targetDensityFunction) != len(subdivRegionInfo):
        raise GridException("When specifying a list of functions, as many as subdivision regions should be specified")

    if not excludeFunction:
        excludeFunction = _defaultExcludeFunction
    if not checkSubDivFunction:
        checkSubDivFunction = _defaultSubDivFunction

    if minSquares >= maxSquares:
        raise GridException("Minimal number of cells must be smaller than maximum")
    
    if excludeFunction is None and maxSquares < sum([ x["startsubdiv"]**2 for x in subdivRegionInfo ]):
        raise GridException("Requested initial subdivision leads to more cells than the specified maximum")

    # Note that this isn't really the total mass, we should multiply with Dd**2 for this
    minPointDists = [ x["size"]/float(maxIntegrationSubDiv) for x in subdivRegionInfo ]
    totalMasses = [ _integrateGridFunction(f, x["center"], x["size"], minPointDist) for f, x, minPointDist in zip(targetDensityFunction, subdivRegionInfo, minPointDists) ]
    totalMass = sum(totalMasses)
    #print("total mass in area:", totalMass*lensInfo["Dd"]**2)

    subDivFraction = 0.2
    diffFrac = 2.0
    maxDiffFracIt = 50
    maxSubDivIt = 1000

    massCache = { }

    totalCells = lambda allGrids: sum([ len(grid["cells"]) for grid in allGrids ])

    for i in range(maxDiffFracIt):
        diffFrac /= 2.0

        maxIt = maxSubDivIt
        allGrids = None

        while not allGrids or (totalCells(allGrids) < minSquares and maxIt > 0):
            maxIt -= 1
            subDivFraction /= (1.0 + diffFrac)
            allGrids = _createSubdivisionGridForThreshold_Multi(targetDensityFunction, subdivRegionInfo, totalMass*subDivFraction,
                                                      excludeFunction, maxIntegrationSubDiv, keepLarger, massCache,
                                                      checkSubDivFunction, maxSquares)

        if totalCells(allGrids) < minSquares:
            raise GridException("Unable to find a grid that has a larger number of squares than the specified minimum")

        if totalCells(allGrids) <= maxSquares:
            return allGrids

        diffFrac /= 2.0
        maxIt = maxSubDivIt

        while totalCells(allGrids) > maxSquares and maxIt > 0:
            maxIt -= 1
            subDivFraction *= (1.0 + diffFrac)
            allGrids= _createSubdivisionGridForThreshold_Multi(targetDensityFunction, subdivRegionInfo, totalMass*subDivFraction,
                                                      excludeFunction, maxIntegrationSubDiv, keepLarger, massCache,
                                                      checkSubDivFunction, maxSquares)

        if totalCells(allGrids) > maxSquares:
            raise GridException("Unable to find a grid that has a smaller number of squares than the specified maximum")

        if totalCells(allGrids) >= minSquares:
            return allGrids

    raise GridException("A grid of the requested size could not be constructed")

def createSubdivisionGridForFunction(targetDensityFunction, size, center, minSquares, maxSquares, startSubDiv = 1, 
        excludeFunction = None,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        checkSubDivFunction = None):
    """Creates a grid of which the cell density is based on the values provided by a function `f`,
    and do the refinement in such a way that the resulting number of grid cells is within the
    specified bounds. This is the core function that's used by :func:`createSubdivisionGrid`
    and :func:`createSubdivisionGridForFITS`, which themselves may me more straightforward to 
    use for specific purposes.
    
    Arguments:

     - `targetDensityFunction`: a function which will be called with a (Ny,Nx,2) shaped NumPy array, describing
       X,Y coordinates on a 2D grid. The function must return an (Ny, Nx) shaped NumPy array
       containing the densities at the specified coordinates. Note that because the algorithm
       below will numerically integrate cells in the grid, it is probably not a good idea to
       let this function return negative values.
     - `size` and `center`: specifies the center X,Y coordinates of the grid, as well as the
       grid's width and height.
     - `minSquares` and `maxSquares`: stop looking for the required grid refinement when the
       total number of grid cells lies between these values.
     - `startSubDiv`: by default, a single, large grid cell will be subdivided until the required
       refinement has been reached. By setting this option to 5 for example, you'd start from
       a uniform 5x5 grid, of which each cell will be refined further.
     - `excludeFunction`: if specified, for each grid cell it will be called as 
       ``excludeFunction(cellCenter, cellSize, True)`` for the initial uniform grid (see `startSubDiv`), and
       as ``excludeFunction(cellCenter, cellSize, False)`` for the next subdivision steps. If the function returns ``True``
       then that particular cell will not be included in the final grid.
     - `maxIntegrationSubDiv`: to determine if a cell needs to be subdivided, the algorithm
       will integrate the contents numerically. This number specifies the maximum number of
       parts that are used for the entire `size` range. Smaller cells will be subdivided in
       fewer parts according to their size.
     - `keepLarger`: if ``True``, after subdividing a grid cell, the original cell will also
       be kept in the list of cells.
     - `checkSubDivFunction`: TODO

    Refinement algorithm: in a helper routine, each cell will be split into four smaller cells
    whenever the integrated value of the cell exceeds some threshold. The main function then
    iteratively looks for an appropriate threshold, so that the final number of cells lies
    within the desired bounds.
    """

    return createMultiSubdivisionGridForFunction(targetDensityFunction,
                [{"size": size, "center": center, "startsubdiv": startSubDiv}],
                minSquares, maxSquares, excludeFunction, maxIntegrationSubDiv, keepLarger,
                checkSubDivFunction)[0]

def createMultiSubdivisionGrid(subdivRegionInfo, # [ { size, center, startsubdiv } ]
        lensInfo, minSquares, maxSquares, 
        excludeFunction = None,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        ignoreOffset = True,
        useAbsoluteValues = True,
        checkSubDivFunction = None):
    """TODO"""
    
    if not excludeFunction:
        excludeFunction = _defaultExcludeFunction
    if not checkSubDivFunction:
        checkSubDivFunction = _defaultSubDivFunction

    if type(lensInfo) != list:
        allLensInfo = [ lensInfo for x in subdivRegionInfo ]
    else:
        if len(lensInfo) != len(subdivRegionInfo):
            raise GridException("Expecting the same amount of lensInfo instances as regions")
        allLensInfo = lensInfo
        
    allFunctions = []
    for lensInfo in allLensInfo:
        densPoints = lensInfo.getDensityPoints()
        if useAbsoluteValues:
            densPoints = np.absolute(densPoints)

        offset = densPoints.min()
        if offset < 0:
            raise GridException("Detected negative value in density points, can't handle this")

        if ignoreOffset:
            densPoints = densPoints - offset

        f = gridfunction.GridFunction(densPoints, lensInfo.getBottomLeft(), lensInfo.getTopRight())
        allFunctions.append(f.evaluate)

    return createMultiSubdivisionGridForFunction(allFunctions, subdivRegionInfo, minSquares, maxSquares,
                                        excludeFunction, maxIntegrationSubDiv, keepLarger, checkSubDivFunction)

def createSubdivisionGrid(size, center, lensInfo, minSquares, maxSquares, startSubDiv = 1, 
        excludeFunction = None,
        maxIntegrationSubDiv = 256,
        keepLarger = False,
        ignoreOffset = True,
        useAbsoluteValues = True,
        checkSubDivFunction = None):

    """Creates a subdivision grid that's based on the specified gravitational lens information.
    Internally, the function :func:`createSubdivisionGridForFunction` is called, where you 
    can find the explanation of the algorithm used.

    Arguments:

     - `size`: the width and height of the grid, will be used as the `size` parameter
       in :func:`createSubdivisionGridForFunction`
     - `center`: the center of the grid in the re-centered coordinates. Will be used
       as the `center` parameter in :func:`createSubdivisionGridForFunction`.
     - `lensInfo`: information about the gravitational lens, see
       :class:`LensInfo<grale.plotutil.LensInfo>` (or
       :class:`DensInfo<grale.plotutil.DensInfo>`).
     - `minSquares`, `maxSquares`, `startSubDiv`, `excludeFunction`, `maxIntegrationSubDiv`,
       `keepLarger`: see :func:`createSubdivisionGridForFunction`
     - `ignoreOffset`: adds or subtracts a value from the data so that the minimum value
       becomes zero.
     - `useAbsoluteValues`: set to ``False`` to allow negative values (probably not a good
       idea).
     - `checkSubDivFunction`: TODO
    """
    return createMultiSubdivisionGrid([{ "size": size, "center": center, "startsubdiv": startSubDiv}],
                                      lensInfo, minSquares, maxSquares, excludeFunction, maxIntegrationSubDiv,
                                      keepLarger, ignoreOffset, useAbsoluteValues,
                                      checkSubDivFunction)[0]

def fitMultiplePlummerLens(gridCells, Dd, targetDensityFunction, sizeFactor = 1.7):
    """Fits a mass distribution based on multiple Plummer base functions that are
    arranged according to specified grid cells, to a specific target function. For
    each cell, the center of the cell and the centers of each side will be used to
    evaluate the target function. A least squares approach is then used to find the
    weigths of the basis functions that minimize the differences at these points.
    
    The function returns an instance of the :class:`MultiplePlummerLens<grale.lenses.MultiplePlummerLens>`
    class.

    Examples can be found in the notebook `fittest.ipynb <_static/fittest.ipynb>`_.

    Arguments:

     - `gridCells`: the grid to base the layout of the Plummer basis functions on
     - `Dd`: the angular diameter distance to be used when creating the
       :class:`MultiplePlummerLens<grale.lenses.MultiplePlummerLens>` instance.
     - `targetDensityFunction`: a function which will be called with a (Ny,Nx,2) shaped NumPy array, describing
       X,Y coordinates on a 2D grid. The function must return an (Ny, Nx) shaped NumPy array
       containing the densities at the specified coordinates. Note that because the algorithm
       below will numerically integrate cells in the grid, it is probably not a good idea to
       let this function return negative values.
     - `sizeFactor`: a factor to be used to convert the size of a grid cell to
       the width of the Plummer basis function. For a cell with a size of 1 arcsec,
       a Plummer width of 1.7 arcsec will be used.
    """

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

class MultiGridCreator(object):
    """TODO"""

    def __init__(self, regionSize = None, regionCenter = None, multiRegionInfo = None,
                 defaultRandomFraction = 0.05):
        """TODO""" # defaultRandomFraction can also be a function
      
        self.defaultRandomFraction = defaultRandomFraction
        self.multiRegionInfo = []
        self.setRegionSize(regionSize, regionCenter, multiRegionInfo)

    def setRegionSize(self, regionSize = None, regionCenter = None, 
                      multiRegionInfo = None):
        """TODO"""

        if multiRegionInfo:
            if not (regionSize is None and regionCenter is None):
                raise GridException("If multiRegionInfo is set, regionSize and regionCenter must be None")

            self.multiRegionInfo = []
            if type(multiRegionInfo) == dict:
                multiRegionInfo = [ multiRegionInfo ]

            for reg in multiRegionInfo:
                ctr = [0.0,0.0] if (not "center" in reg or not reg["center"]) else copy.copy(reg["center"])
                sz = float(reg["size"])
                if sz <= 0:
                    raise GridException("Invalid region size {} was specified".format(sz))

                self.multiRegionInfo.append({"center": ctr, "size": sz})

        else: # multiRegionInfo is None
            if not regionSize:
                raise GridException("No region info was specified")

            ctr = [0.0, 0.0] if not regionCenter else copy.copy(regionCenter)
            sz = float(regionSize)
            if sz <= 0:
                raise GridException("Invalid region size {} was specified".format(sz))

            self.multiRegionInfo = [ { "size": sz, "center": ctr }]

    def _getRandomizedRegions(self, randomFraction):
        randomizedMultiReg = []
        for reg in self.multiRegionInfo:

            if type(randomFraction) == float: # Just a number, use the default method
                w = reg["size"]
                dx = (random.random()-0.5) * w*randomFraction
                dy = (random.random()-0.5) * w*randomFraction
                c = [ reg["center"][0] + dx, reg["center"][1] + dy ]
            else: # assume it's something we can call to obtain the grid dimensions
                w, c = randomFraction(reg["size"], copy.copy(reg["center"]))

            randomizedMultiReg.append({"size": w, "center": c})

        return randomizedMultiReg

    def getUniformGrid(self, subDiv, randomFraction = None, excludeFunction = None):
        """TODO"""
        
        if randomFraction is None:
            randomFraction = self.defaultRandomFraction

        randomizedMultiReg = self._getRandomizedRegions(randomFraction)
        if type(subDiv) != list and type(subDiv) != tuple:
            subDiv = [ subDiv for x in randomizedMultiReg ] # One subdiv for each multi-grid part

        if len(randomizedMultiReg) != len(subDiv):
            raise Exception("There are {} multi-grid entries, but you specified {} subdivisions".format(len(randomizedMultiReg),len(subDiv)))

        for part,sd in zip(randomizedMultiReg,subDiv):
            part["axissubdiv"] = sd
        return createMultiUniformGrid(randomizedMultiReg, excludeFunction)

    def getSubdivisionGrid(self, lensOrLensInfo, minSquares, maxSquares, startSubDiv = 1, 
                           randomFraction = None, lensFilter = None, lensFilterParams = {},
                           lensInfoFilter = None, lensInfoFilterParams = {},
                           excludeFunction=None, checkSubDivFunction=None,
                           renderer = "default", feedbackObject = "default"):
        """TODO"""
        
        if lensFilter is None:
            lensFilter = lambda x, bl, tr, **params: x
        if lensInfoFilter is None:
            lensInfoFilter = lambda x, **params: x

        randomizedMultiReg = self._getRandomizedRegions(randomFraction)

        if type(startSubDiv) == list or type(startSubDiv) == tuple: # A list with a subdiv entry for each region
            if len(startSubDiv) != len(randomizedMultiReg):
                raise GridException("Expecting as many start subdivisions as regions")
        else:
            # Just a number, create a list 
            startSubDiv = [ startSubDiv for x in randomizedMultiReg ]

        for part, sd in zip(randomizedMultiReg, startSubDiv):
            part["startsubdiv"] = sd

        from . import plotutil

        if isinstance(lensOrLensInfo, plotutil.DensInfo):
            lensInfo = lensOrLensInfo
            lensInfo = lensInfoFilter(lensInfo, **lensInfoFilterParams)
        else:
            # Treat each grid part separately
            allLensInfos = []
            for part in randomizedMultiReg:
                w,c = part["size"], part["center"]

                extraFrac = 1.0001
                
                bl = [ c[0] - (w*extraFrac)/2, c[1] - (w*extraFrac)/2 ]
                tr = [ c[0] + (w*extraFrac)/2, c[1] + (w*extraFrac)/2 ]
                lensOrLensInfo = lensFilter(lensOrLensInfo, bl, tr, **lensFilterParams)

                lensInfo = plotutil.LensInfo(lensOrLensInfo, bottomleft=bl, topright=tr)
                lensInfo.getDensityPoints(renderer, feedbackObject)
            
                lensInfo = lensInfoFilter(lensInfo, **lensInfoFilterParams)
                allLensInfos.append(lensInfo)

            lensInfo = allLensInfos

        return createMultiSubdivisionGrid(randomizedMultiReg, lensInfo, minSquares, maxSquares,
                                                     excludeFunction=excludeFunction,
                                                     checkSubDivFunction=checkSubDivFunction)

def main():
    #grid = createUniformGrid(123, [-1, -2], 5)
    #pprint.pprint(grid)
    #g = _fractionalGridToRealGrid(grid)
    #_debugPlot(g)

    import grale.lenses as lenses
    import grale.plotutil as plotutil
    from grale.constants import ANGLE_ARCSEC
    import matplotlib.pyplot as plt


    U = ANGLE_ARCSEC
    lensInfo = plotutil.LensInfo(**{
        "lens": lenses.GravitationalLens.load("/home/jori/projects/grale-lensmodels-git/2012MNRAS.425.1772L/reallens.lensdata"),
        "bottomleft": [-100*U, -100*U],
        "topright": [100*U, 100*U],
        })
    plotutil.plotDensity(lensInfo, feedbackObject = "none")
    plt.show()
    #pprint.pprint(lensInfo)

    grid = createSubdivisionGrid(190*U, [0,0], lensInfo, 900, 1000, startSubDiv=7)
    print("# num cells", len(grid["cells"]))
    #pprint.pprint(grid)

    g = _fractionalGridToRealGrid(grid)
    _debugPlot(g)
    plotutil.plotSubdivisionGrid(grid)
    plt.show()

if __name__ == "__main__":
    main()
