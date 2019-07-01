import numpy as np

def _findBestPermutation(observed, predictions):

    def findClosest(point, arr):
        bestDist, bestIdx = float("inf"), None

        for i in range(len(arr)):
            dist2 = np.sum((arr[i]-point)**2)

            if dist2 < bestDist:
                bestDist, bestIdx = dist2, i

        if bestIdx is None:
            raise Exception("Couldn't find closest point!")

        return bestDist**0.5, bestIdx

    bestDist, bestPerm = float("inf"), None

    positions = [ { "pos": x} for x in observed ]
    for i in range(len(positions)):
        positions[i]["origidx"] = i

    import itertools

    for posPerm in itertools.permutations(positions):
        predCopy = predictions[:] # we make a copy because we're going to remove things
        predPerm = [ None ] * len(positions)

        curDist = 0.0
        for point in posPerm:
            pos = point["pos"]

            dist, idx = findClosest(pos, predCopy)
            predPerm[point["origidx"]] = predCopy[idx]
            del predCopy[idx]

            curDist += dist

            if curDist >= bestDist: # can't get any better
                break

        if curDist < bestDist:
            bestDist, bestPerm = curDist, predPerm

            for x in bestPerm:
                if x is None:
                    raise Exception("UNEXPECTED!")

    return bestDist, bestPerm

def _align(pred, obs, maxPermSize):

    if len(pred) < len(obs):
        raise Exception("Not enough predictions to match observations")

    pred = [ np.array(p) for p in pred ]
    obs = [ np.array(o) for o in obs ]

    distances = [ [] for i in range(len(obs))]
    for i in range(len(obs)):
        theta_obs = obs[i]

        for j in range(len(pred)):
            theta_pred = pred[j]
            distSquared = np.sum((theta_pred-theta_obs)**2)
            distances[i].append((distSquared,j))

    havePoints = set()
    needDetailedComparision = False
    # Sort each distance list to reveal the closest point
    result = [ None for o in obs ]

    for i in range(len(distances)):
        l = distances[i]
        l.sort()
        closestPredPoint = l[0][1]
        if closestPredPoint in havePoints:
            needDetailedComparision = True
            break
        else:
            havePoints.add(closestPredPoint)
            result[i] = pred[closestPredPoint]

    if not needDetailedComparision:
        return result

    if len(obs) > maxPermSize:
        raise Exception(f"Need more detailed comparision, but number of observations {len(obs)} exceeds {maxPermSize}; you can increase 'maxPermSize', but it may take some time.")

    # Go over permurations
    bestDist, bestPerm = _findBestPermutation(obs, pred)
    return bestPerm

def calculateImagePredictions(imgList, lensModel, cosmology=None,
                     reduceImages="average",
                     useAverageBeta=True,
                     maxPermSize=7):

    from . import multiplane
    from . import privutil
    from numpy.linalg import inv

    if hasattr(lensModel, "getLensDistance"):
        zd = min(cosmology.findRedshiftForAngularDiameterDistance(lensModel.getLensDistance()))
        lensModel = [ (lensModel, zd) ]

    useTrace = True
    if type(lensModel) == multiplane.MultiLensPlane:
        lensPlane = lensModel
    elif type(lensModel) == multiplane.MultiImagePlane:
        lensPlane = lensModel.getLensPlane()
    else:
        lensPlane = multiplane.MultiLensPlane(lensModel, [0,0], [0,0], 1, 1, None, "none", cosmology)
        useTrace = False

    if useTrace:
        if cosmology:
            raise Exception("When using a multi-lensplane to trace positions, the internal cosmological model is used; this parameter must be set to None")
    else:
        cosmology = privutil.initCosmology(cosmology)
        if not cosmology:
            raise Exception("Cosmological model is not set")

    def getImagePoint_Average(imgDat, imgNum):
        pts = np.array([imgDat.getImagePointPosition(imgNum, p) for p in range(imgDat.getNumberOfImagePoints(imgNum))])
        return np.average(pts,0)

    getImagePoint = getImagePoint_Average if reduceImages == "average" else reduceImages

    allPoints = [ ]
    for i in imgList:
        imgDat = i["imgdata"]
        z = i["z"]

        imgPlane = multiplane.MultiImagePlane(lensPlane, z)
        imgPos = [ ]
        for imgNum in range(imgDat.getNumberOfImages()):
            theta = getImagePoint(imgDat, imgNum)

            if useTrace:
                beta = imgPlane.traceTheta(theta)
                imgPos.append({ "theta": theta, "beta": beta })
            else:
                beta, betaDerivs = imgPlane.getBetaAndDerivatives(theta)
                imgPos.append({ "theta": theta, "beta": beta,
                                "betaderivs": betaDerivs, "invbetaderivs": inv(betaDerivs)})

        allPoints.append((imgPos, z))

    sources = []
    for imgPos, z in allPoints:

        sourceInfo = [ ]
        betas = [ np.average([d["beta"] for d in imgPos], 0) ] if useAverageBeta else [ d["beta"] for d in imgPos ]

        if useTrace:
            imgPlane = multiplane.MultiImagePlane(lensPlane, z)

            observedThetas = [ x["theta"] for x in imgPos ]
            predictedThetas = [ [] for i in range(len(observedThetas)) ]
            # Trace each beta, and align the resulting position with the observations
            for beta in betas:
                thetas = imgPlane.traceBetaApproximately(beta)
                thetas = _align(thetas, observedThetas, maxPermSize)

                for i in range(len(thetas)):
                    predictedThetas[i].append(thetas[i])

            for i in range(len(observedThetas)):
                sourceInfo.append({ "theta_obs": observedThetas[i],
                                    "theta_pred": predictedThetas[i] })
        else:
            for x in imgPos:
                thetaPred = [ ]
                theta0, beta0, invderiv0 = x["theta"], x["beta"], x["invbetaderivs"]
                for beta1 in betas:
                    dbeta = beta1-beta0
                    dtheta = invderiv0.dot(dbeta)
                    thetaPred.append(theta0 + dtheta)

                sourceInfo.append({ "theta_obs": theta0,
                                    "theta_pred": thetaPred })

        sources.append(sourceInfo)
    return sources

def _getPredictions(theta_pred, avgImage):
    if not avgImage:
        return theta_pred
    return [ np.average([t for t in theta_pred], 0) ]

def calculateRMS_all(predictions, avgImage = False):
    sum2 = 0
    count = 0
    for sourceInfo in predictions:
        for imgInfo in sourceInfo:
            theta_obs = imgInfo["theta_obs"]

            for theta_pred in _getPredictions(imgInfo["theta_pred"], avgImage):
                dt = theta_pred-theta_obs
                sum2 += dt[0]**2 + dt[1]**2
                count += 1

    return (sum2/count)**0.5

def calculateRMS_averageRMSPerSource(predictions, avgImage = False):

    sourceSum = 0

    for sourceInfo in predictions:
        sum2 = 0
        count = 0
        for imgInfo in sourceInfo:
            theta_obs = imgInfo["theta_obs"]

            for theta_pred in _getPredictions(imgInfo["theta_pred"], avgImage):
                dt = theta_pred-theta_obs
                sum2 += dt[0]**2 + dt[1]**2
                count += 1

        sourceSum += (sum2/count)**0.5

    return sourceSum/len(predictions)

def calculateRMS_averageSquaredPerSource(predictions, avgImage = False):

    sourceSum = 0

    for sourceInfo in predictions:
        sum2 = 0
        count = 0
        for imgInfo in sourceInfo:
            theta_obs = imgInfo["theta_obs"]
            for theta_pred in getPredictions(imgInfo["theta_pred"], avgImage):
                dt = theta_pred-theta_obs
                sum2 += dt[0]**2 + dt[1]**2
                count += 1

        sourceSum += sum2/count

    return (sourceSum/len(predictions))**0.5

