
# Quite straightforward port of a C++ routine
def createTriangulationFromHull(pts, maxImgPoints = 50, maxRelativeDistance = 0.02):

    def getScale():
        X, Y = pts[:,0], pts[:,1]
        dx, dy = X.max()-X.min(), Y.max()-Y.min()
        return (dx**2+dy**2)**0.5

    def getDistance(p1, p2, p0):
        x0, y0 = p0
        x1, y1 = p1
        x2, y2 = p2

        dx, dy = x2-x1, y2-y1
        denom = (dx**2+dy**2)**0.5
        return abs(dy*x0-dx*y0+x2*y1-y2*x1)/denom

    idxOneThird = len(pts)//3
    idxTwoThirds = len(pts)*2//3

    if not (0 < idxOneThird and idxOneThird < idxTwoThirds and idxTwoThirds < len(pts)):
        raise Exception("Couldn't find three separate points to start from")

    imgPoints = [ 0, idxOneThird, idxTwoThirds, len(pts) ]

    scale = getScale()
    maxToleratedDistance = scale*maxRelativeDistance

    while len(imgPoints) < maxImgPoints:

        allDistOk, newImgPoints = True, [ 0 ]

        for i in range(len(imgPoints) - 1):

            if len(newImgPoints) >= maxImgPoints:
                break

            idx0, idx1 = imgPoints[i:i+2]
            if idx0 >= idx1:
                raise Exception("Internal error: idx0 >= idx1")

            if idx1 > idx0+1:
                p0, p1 = pts[idx0,:], pts[idx1%len(pts),:]

                maxDist, maxDistIdx = 0, -1
                for idx in range(idx0, idx1):
                    testPoint = pts[idx,:]
                    dist = getDistance(p0, p1, testPoint)
                    if dist > maxDist:
                        maxDist, maxDistIdx = dist, idx

                if maxDist > maxToleratedDistance:
                    newImgPoints.append(maxDistIdx)
                    allDistOk = False

            newImgPoints.append(idx1)

        imgPoints = newImgPoints
        if allDistOk:
            break

    imgPoints = imgPoints[:-1]

    from scipy.spatial import Delaunay
    tri = Delaunay(pts[imgPoints])

    from shapely.geometry import Point
    from shapely.geometry.polygon import Polygon

    simp = [ ]

    originalPoly = Polygon(pts)
    for p0, p1, p2 in tri.simplices:
        xc, yc = (pts[imgPoints[p0]] + pts[imgPoints[p1]] + pts[imgPoints[p2]])/3.0
        
        if originalPoly.contains(Point(xc, yc)):
            simp.append([p0, p1, p2])

    return imgPoints, simp
