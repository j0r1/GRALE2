"""Module meant for various utilities, but currently only RMS calculation."""

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

def _findRealTheta(imgPlane, beta, theta):
    from .constants import ANGLE_ARCSEC
    import scipy.optimize as opt

    startDiff = (imgPlane.traceTheta(theta) - beta)/ANGLE_ARCSEC

    def f(x):
        return (imgPlane.traceTheta(x) - beta)/ANGLE_ARCSEC

    r = opt.fsolve(f, x0=theta)
    #endDiff = (imgPlane.traceTheta(r) - beta)/ANGLE_ARCSEC
    #print("start diff", sum(startDiff**2)**0.5, "arcsec")
    #print("end diff", sum(endDiff**2)**0.5, "arcsec")
    return r

def calculateImagePredictions(imgList, lensModel, cosmology=None,
                     reduceImages="average",
                     useAverageBeta=True,
                     maxPermSize=7,
                     useFSolve=True):
    r"""For a set of images and a lens model, the predicted images for
    each observed image are calculated. The results can subsequently
    be used in the :func:`calculateRMS` function to obtain the RMS.

    Parameters:
     - `imgList`: a list of dictionaries with entries

         - ``"imgdata"``: :class:`ImagesData <grale.images.ImagesData>` instance
           that describes the observed images.
         - ``"z"``: the redshift for these observed images.

     - `lensModel`: this can be a :class:`MultiLensPlane <grale.multiplane.MultiLensPlane>`
       instance, a :class:`MultiImagePlane <grale.multiplane.MultiImagePlane>` instance,
       a list of (:class:`lens model <grale.lenses.GravitationalLens>`, redshift)
       tuples, or just a single :class:`model<grale.lenses.GravitationalLens>`.

       In case only one or more lens models are specified, the predicted image 
       positions are estimated by calculating the derivatives of :math:`\vec{\beta}(\vec{\theta})`
       and using them to compensate for small differences in source plane
       positions :math:`\vec{\beta}`; or if `useFSolve` is ``True``, several
       iterations will be used to approximate the true solution.

       On the other hand, if e.g. a MultiLensPlane
       instance is used, then the :func:`traceBetaApproximately <grale.multiplane.MultiImagePlane.traceBetaApproximately>`
       function is called to obtain estimates of image plane positions corresponding
       to a source plane position. This is more computationally demanding, but should
       yield a more correct result.

     - `cosmology`: must only be specified somehow (e.g. "default" if a default cosmology
       was set if a (list of) lens models was specified. If e.g. a MultiLensPlane was
       used, the internally stored cosmological model will be used instead.

     - `reduceimages`: each image in an images data set will be converted to a single
       point. By default, this is done by averaging the image point positions. If something
       else needs to be done, you can specify a function here, which will be called with
       two arguments: the :class:`ImagesData <grale.images.ImagesData>` instance and the
       image index corresponding to the particular image.

     - `useAverageBeta`: if ``True``, then the back-projected image points are averaged
       to estimate a single source position. If ``False``, then each back-projected image
       is used as an estimate of the source position.

     - `maxPermSize`: when a source plane position is used to estimate the corresponding
       image plane positions, these predicted positions are grouped with the observed
       positions. If the derivatives-based method is used, then this is automatically
       possible, but if the more accurate tracing is used, we need to find out which
       predicted point corresponds to which observed point. As a first attempt, the
       procedure tries to use the points closest to the observed ones, but this may not
       be possible and different permutations of the positions will then be explored.
       Since this can become quite slow, an error is generated if more elements than
       this would need to be permutated.

     - `useFSolve`: if only a lens model is used, this cause SciPy's `fsolve <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.fsolve.html>`_
       to be used to zoom in on the real image plane vectors. Otherwise, a one-step
       approach is used based on the deflection angle derivatives.
    """

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
                imgPos.append({ "theta": theta, "beta": beta, "z": z,
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
                                    "theta_pred": predictedThetas[i],
                                    "beta_est": betas})
        else:
            for x in imgPos:
                if useFSolve:
                    imgPlane = multiplane.MultiImagePlane(lensPlane, x["z"])

                thetaPred = [ ]
                theta0, beta0, invderiv0 = x["theta"], x["beta"], x["invbetaderivs"]
                for beta1 in betas:
                    if useFSolve:
                        # We're looking for the images of beta1, starting from the
                        # estimate at observed theta0
                        thetaPred.append(_findRealTheta(imgPlane, beta1, theta0)) 
                    else:
                        dbeta = beta1-beta0
                        dtheta = invderiv0.dot(dbeta)
                        thetaPred.append(theta0 + dtheta)

                sourceInfo.append({ "theta_obs": theta0,
                                    "theta_pred": thetaPred,
                                    "beta_est": betas })

        sources.append(sourceInfo)
    return sources

def _getPredictions(theta_pred, avgImage):
    if not avgImage:
        return theta_pred
    return [ np.average([t for t in theta_pred], 0) ]

def calculateRMS(predictions, angularUnit, avgImage = False):
    r"""Using the output of :func:`calculateImagePredictions`, the RMS value
    is calculated. 
    Using the notation from Appendix A of `Williams et al (2018) <https://ui.adsabs.harvard.edu/abs/2018MNRAS.480.3140W/abstract>`_,
    in general there are :math:`i=1..I` sources, with :math:`j=1..J_i` images each. When using
    each back-projected observed image position :math:`\vec{\theta}_{i,j}`
    as a separate estimate of the source position, there will be :math:`k=1..J_i`
    image predictions :math:`\vec{\theta}_{i,j,k}` for each observed position.

    When treating each predicted-observed difference of all sources equally,
    the RMS will be given by

    .. math::

        {\rm RMS}_{\rm full}^2 = \frac{1}{K}\sum_{i=1}^I \sum_{j=1}^{J_i} \sum_{k=1}^{J_i} \left|\vec{\theta}_{i,j,k} - \vec{\theta}_{i,j}\right|^2

    where :math:`K = \sum_{i=1}^{I} J_i^2` is the total number of terms in this summation.
    Alternatively, we could take the average of the right-most sum, thereby averaging
    the squared differences for each observed image point. Sources with more images
    will still have more terms in the rest of the summation, but only one per image
    (and not :math:`J_i` for each image). This results in

    .. math::

        {\rm RMS}_{\rm equalimages}^2 = \frac{1}{J}\sum_{i=1}^I \sum_{j=1}^{J_i} \frac{1}{J_i} \sum_{k=1}^{J_i} \left|\vec{\theta}_{i,j,k} - \vec{\theta}_{i,j}\right|^2

    Here, :math:`J = \sum_{i=1}^I J_i`. We could also average all :math:`J_i^2` terms
    per source first, and then average over the remaining :math:`I` terms, one for each
    source:

    .. math::

        {\rm RMS}_{\rm equalsources}^2 = \frac{1}{I}\sum_{i=1}^I \frac{1}{J_i^2} \sum_{j=1}^{J_i} \sum_{k=1}^{J_i} \left|\vec{\theta}_{i,j,k} - \vec{\theta}_{i,j}\right|^2

    This function returns a dictionary with these three RMS values (non-squared),
    expressed in the specified `angularUnit`. 
    
    If `avgImage` is ``True``, then all image predictions for each observed position 
    are averaged out first. This will only make sense in case `useAverageBeta` was ``False``
    in the call to :func:`calculateImagePredictions` since otherwise there's only
    one predicted position.

    Depending on the settings of such flags, and depending on the use of the
    derivatives-based prediction or the more accurate tracing, some of these results
    can yield the same value.
    """

    fullRms, equalSourceRMS, equalImageRMS = [], [], []

    for sourceInfo in predictions:

        sumImages = [ ]

        for imgInfo in sourceInfo:
            theta_obs = imgInfo["theta_obs"]
            sumPredictedPerObserved =  [ ]

            predThetaList = _getPredictions(imgInfo["theta_pred"], avgImage)
            for theta_pred in predThetaList:
                dt = theta_pred-theta_obs
                dist2 = dt[0]**2 + dt[1]**2

                sumPredictedPerObserved.append(dist2)
                sumImages.append(dist2)
                fullRms.append(dist2)

            equalImageRMS.append(np.average(sumPredictedPerObserved))

        equalSourceRMS.append(np.average(sumImages))

    return {
        "full": np.average(fullRms)**0.5/angularUnit,
        "equalsources": np.average(equalSourceRMS)**0.5/angularUnit,
        "equalimages": np.average(equalImageRMS)**0.5/angularUnit
    }

