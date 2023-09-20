"""Module meant for various utilities, but currently only RMS calculation."""

import numpy as np
from .inverters import _getNumHelpers
from multiprocessing import Pool
import copy

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

    #startDiff = (imgPlane.traceTheta(theta) - beta)/ANGLE_ARCSEC

    def f(x):
        return (imgPlane.traceTheta(x) - beta)/ANGLE_ARCSEC

    r = opt.fsolve(f, x0=theta)
    #endDiff = (imgPlane.traceTheta(r) - beta)/ANGLE_ARCSEC
    #print("start diff", sum(startDiff**2)**0.5, "arcsec")
    #print("end diff", sum(endDiff**2)**0.5, "arcsec")
    return r

class _ImagePlaneWrapper(object):
    def __init__(self, imgPlane, lensPlane):
        self.imgPlane = imgPlane
        self.lensPlane = lensPlane
        
    def traceTheta(self, theta):
        lens = self.lensPlane.getLens()
        Dds = self.imgPlane.getDds()
        Ds = self.imgPlane.getDs()
        return lens.traceTheta(Ds, Dds, theta)
        
    def getBetaAndDerivatives(self, theta):
        lens = self.imgPlane.getLens()
        Dds = self.imgPlane.getDds()
        Ds = self.imgPlane.getDs()
        beta = lens.traceTheta(Ds, Dds, theta)
        axx, ayy, axy = lens.getAlphaVectorDerivatives(theta) * (Dds/Ds)
        derivs = np.array([[1-axx, -axy],[-axy, 1-ayy]], dtype=np.double)
        return beta, derivs
        
    def traceBetaApproximately(self, beta):
        return self.imgPlane.traceBeta(beta)

def _commonInitFindOptRetrace(imgList, lensModel, cosmology, reduceImages):

    from . import multiplane
    from . import privutil
    from . import images
    from numpy.linalg import inv

    cosmology = privutil.initCosmology(cosmology)
    origLensModel = lensModel

    def findZ(lensModel):
        return min(cosmology.findRedshiftForAngularDiameterDistance(lensModel.getLensDistance()))

    if hasattr(lensModel, "getLensDistance"):
        if not cosmology:
            raise Exception("When using a lens model, the cosmological model should be set")
        zd = findZ(lensModel)
        lensModel = [ (lensModel, zd) ]

    useTrace = True
    if type(lensModel) == multiplane.MultiLensPlane:
        lensPlane = lensModel
        createImgPlaneFn = lambda lp, zs: multiplane.MultiImagePlane(lp, zs)
    elif type(lensModel) == multiplane.MultiImagePlane:
        lensPlane = lensModel.getLensPlane()
        createImgPlaneFn = lambda lp, zs: multiplane.MultiImagePlane(lp, zs)
    elif type(lensModel) == images.LensPlane:
        if not cosmology:
            raise Exception("When using a regular lens plane to trace positions, the cosmological model should be set")

        lensPlane = lensModel
        zd = findZ(lensPlane.getLens())
        createImgPlaneFn = lambda lp, zs: _ImagePlaneWrapper(images.ImagePlane(lensPlane, cosmology.getAngularDiameterDistance(zs), cosmology.getAngularDiameterDistance(zd, zs)), lensPlane)
    else:
        lensPlane = multiplane.MultiLensPlane(lensModel, [0,0], [0,0], 2, 2, None, "none", cosmology)
        useTrace = False
        createImgPlaneFn = lambda lp, zs: multiplane.MultiImagePlane(lp, zs)

    if useTrace:
        if type(lensModel) == multiplane.MultiLensPlane and cosmology:
            raise Exception("When using a multi-lensplane to trace positions, the internal cosmological model is used; this parameter must be set to None")

    else:
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

        imgPlane = createImgPlaneFn(lensPlane, z)
        imgPos = [ ]
        for imgNum in range(imgDat.getNumberOfImages()):
            theta = getImagePoint(imgDat, imgNum)

            if useTrace:
                beta = imgPlane.traceTheta(theta)
                imgPos.append({ "theta": theta, "beta": beta, "z": z })
            else:
                beta, betaDerivs = imgPlane.getBetaAndDerivatives(theta)
                imgPos.append({ "theta": theta, "beta": beta, "z": z,
                                "betaderivs": betaDerivs, "invbetaderivs": inv(betaDerivs)})

        allPoints.append((imgPos, z))

    return cosmology, lensPlane, origLensModel, useTrace, createImgPlaneFn, allPoints

def nelderMeadSourcePositionOptimizer(imgIdx, imgPos, traceFunction, feedbackObject = None,
                   betaAvgToleranceArcsec = 1e-5, betaAvgPenalty = 10000, nmRounds = 3, nmMaxFev = 50):
    """This is the default optimizer used by :func:`findOptimizedSourcePositions` and
    :func:`parallelFindOptimizedSourcePositions`. It uses the `Nelder-Mead <https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method>`_
    method as `implemented in scipy <https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html>`_.

    Arguments:
     - `imgIdx`, `imgPos`, `traceFunction`: as described in the `optRoutine` argument 
       of :func:`findOptimizedSourcePositions`.
     - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`
     - `betaAvgToleranceArcsec`: depending on the lens model, it is possible that for a
       specific source plane position `beta`, there simply isn't an image plane position
       near the one considered that maps well to this. To filter out such `beta` positions,
       there is this threshold that specifies that if the back-projected (adjusted) image plane
       position still differs more than `betaAvgToleranceArcsec` (so in arc seconds) from
       `beta`, then a penalty should be added.
     - `betaAvgPenalty`: in the case that a penalty (described above) needs to be added,
       so the source plane difference exceeds `betaAvgToleranceArcsec`, this difference
       times `betaAvgPenalty` is added to the function evaluation result.
     - `nmRounds`: the Nelder-Mead routine can be run a number of times, and the best
       result will be kept. This argument specifies that number.
     - `nmMaxFev`: each Nelder-Mead run can have at most this many function evaluations.
    """

    import scipy.optimize as optimize
    from grale.constants import ANGLE_ARCSEC

    feedback = lambda x : None if feedbackObject is None else feedbackObject.onStatus(x)

    if len(imgPos) <= 1:
        raise Exception("There should be at least two image points corresponding to the same source")

    betaStart = np.average([x["beta"] for x in imgPos], axis=0) / ANGLE_ARCSEC
    betaVar = np.sum(( np.std([x["beta"] for x in imgPos], axis=0) / ANGLE_ARCSEC )**2)**0.5

    feedback("Looking for source position for image list index {}".format(imgIdx))

    def f(beta):

        betaPredictionDiffs = []
        thetaPredictionDiffs = []
        for t in [ x["theta"] for x in imgPos ]:
            theta_pred = optimize.fsolve(lambda theta: (traceFunction(theta) - beta*ANGLE_ARCSEC)/ANGLE_ARCSEC, x0=t)
            beta_pred = traceFunction(theta_pred)/ANGLE_ARCSEC

            betaPredictionDiffs.append(np.sum((beta_pred - beta)**2)**0.5)
            thetaPredictionDiffs.append(np.sum(((theta_pred - t)/ANGLE_ARCSEC)**2)**0.5)

        betaPredictionDiffs = np.array(betaPredictionDiffs)
        thetaPredictionDiffs = np.array(thetaPredictionDiffs)

        fitness = 0
        bAvg = np.average(betaPredictionDiffs)
        if bAvg > betaAvgToleranceArcsec:
            fitness += bAvg*betaAvgPenalty

        fitness += np.sum(thetaPredictionDiffs**2)**0.5
        feedback("Fitness for {} is {}".format(beta, fitness))
        return fitness

    if nmRounds < 1:
        raise Exception("Number of Nelder-Mead rounds should be at least 1")

    r = None
    for ridx in range(nmRounds):
        feedback("Starting Nelder-Mead round {}".format(ridx+1))
        betaSimplex = np.array([ betaStart + np.random.normal(0,betaVar/10,size=(2,)) for _ in range(3)])
        r1 = optimize.minimize(f, betaStart, method="Nelder-Mead", options={"initial_simplex": betaSimplex, "maxfev": nmMaxFev})
        if r is None or r1.fun < r.fun:
            r = r1
    feedback("Best source position is {}".format(r.x))
    return r.x * ANGLE_ARCSEC


def findOptimizedSourcePositions(imgList, lensModel, cosmology=None, reduceImages="average", optRoutine=nelderMeadSourcePositionOptimizer, optParams = {}):
    r"""For a given lens model, which typically will not back-project the
    images of a source onto the exact same position, this routine looks for
    a 'best' source position.

    This function returns a list of optimized source positions, one for each 
    entry in `imgList`.

    Arguments:
     - `imgList`: a list of :class:`ImagesData <grale.images.ImagesData>` instances
       describing the images in the strong lensing system. Each image will be reduced
       to a point image position, see also the `reduceImages` argument.
     - `lensModel`: the lens model, which can be represented as a number of things:
    
       - a :class:`GravitationalLens <grale.lenses.GravitationalLens>` instance. In this
         case, the cosmological model `cosmology` must be set as well.
       - a :class: `LensPlane <grale.images.LensPlane>` instance. In this case as well,
         the `cosmology` parameter must be set.
       - a :class:`MultiLensPlane <grale.multiplane.MultiLensPlane>` or
         :class:`MultiImagePlane <grale.multiplane.MultiImagePlane>` instance. The
         `cosmology` parameter must be ``None``, the interally set cosmology will be used.
       - a list of (:class:`GravitationalLens <grale.lenses.GravitationalLens>`, redshift)
         tuples, where the `cosmology` parameter must be set as well.

     - `cosmology`: depending on the way the lens model is specified (see `lensModel`
       argument), this should either be ``None`` or set to the
       :class:`cosmological model <grale.cosmology.Cosmology>` to use.
     - `reduceImages`: each image will need to be reduced to a point image position.
       If set to ``"average"``, then the average position of the image points will be
       used for this. This can also be set to a function, and if so it will be called
       with a :class:`ImagesData <grale.images.ImagesData>` instance as parameter and
       an index of the image. The function should then return a single point that
       corresponds to the specified image.
     - `optRoutine`: this is the optimization routine used to calculate the optimal
       source position for a set of image positions. The default is the
       :func:`nelderMeadSourcePositionOptimizer`. The routine is called with the 
       following arguments: the index into a list of all images (this is useful
       for the parallel version :func:`parallelFindOptimizedSourcePositions`),
       the point image positions for a single source point (a dictionary containing
       ``theta``, ``beta``, ``z`` and depending on the input model also
       ``betaderivs`` and ``invbetaderivs``), the function that's created based on
       the input model to trace an image plane point to its source plane position,
       and the expanded dictionary `optParams` from below.
     - `optParams`: this can be used to pass additional parameters to the
       `optRoutine`.
    """

    if not optRoutine:
        raise Exception("No optimization routine was set")

    cosmology, lensPlane, origLensModel, useTrace, createImgPlaneFn, allPoints = _commonInitFindOptRetrace(imgList, lensModel, cosmology, reduceImages)

    sources = []
    for imgListIdx, (imgPos, z) in enumerate(allPoints):

        imgPlane = createImgPlaneFn(lensPlane, z)
        traceFunction = imgPlane.traceTheta

        beta = optRoutine(imgListIdx, imgPos, traceFunction, **optParams)

        sources.append(beta)

    return sources

def _parallelFindOptimizedSourcePositions_part(part, imgList, lensModel, cosmology, reduceImages, optRoutine, optParams):

    cosmology, lensPlane, origLensModel, useTrace, createImgPlaneFn, allPoints = _commonInitFindOptRetrace(imgList, lensModel, cosmology, reduceImages)

    sources = []
    for imgListIdx, (imgPos, z) in part:
        #t0 = time.time()
        imgPlane = createImgPlaneFn(lensPlane, z)
        traceFunction = imgPlane.traceTheta

        beta = optRoutine(imgListIdx, imgPos, traceFunction, **optParams)
        #t1 = time.time()

        #print("Time for", imgListIdx, "is", t1-t0, "sec, pid = ", os.getpid())
        sources.append([imgListIdx, beta])

    return sources # is only a part


def parallelFindOptimizedSourcePositions(imgList, lensModel, cosmology=None, reduceImages="average", optRoutine=nelderMeadSourcePositionOptimizer, optParams = {}, numThreads=0):
    """This is a version of :func:`findOptimizedSourcePositions` that tries to
    speed up the procedure using multiple threads. It handles several sources in
    parallel.

    Arguments:
     - `imgList`, `lensModel`, `cosmology`, `reduceImages`, `optRoutine` and `optParams` have
       the same meaning as in :func:`findOptimizedSourcePositions`
     - `numThreads`: the number of threads to use to run in parallel. If set to zero or ``None``,
       this will be determined automatically.
    """

    if not optRoutine:
        raise Exception("No optimization routine was set")

    # Check that things work, as well as get allPoints
    cosmology, _, _, _, _, allPoints = _commonInitFindOptRetrace(imgList, lensModel, cosmology, reduceImages)

    numThreads = _getNumHelpers(numThreads)
    pool = Pool(numThreads)
    #print("Created pool for", numThreads, "threads")
    parts = [ [] for i in range(numThreads)]

    for imgListIdx, (imgPos, z) in enumerate(allPoints):
        parts[imgListIdx%numThreads].append([imgListIdx, (imgPos, z)])

    asyncResults = [ pool.apply_async(_parallelFindOptimizedSourcePositions_part, [parts[i], imgList, lensModel, cosmology, reduceImages, optRoutine, optParams]) for i in range(numThreads) ]

    #print("Waiting for everything to finish")
    results = [ r.get() for r in asyncResults ]

    pool.close()
    pool.join()

    sources = [ None for i in range(len(allPoints))]
    for r in results:
        for imgListIdx, beta in r:
            assert sources[imgListIdx] is None, "Unexpected error, writing in same place in output array"
            sources[imgListIdx] = beta

    return sources

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
       tuples, a :class:`LensPlane` <grale.images.LensPlane>` instance,
       or just a single :class:`model<grale.lenses.GravitationalLens>`.

       In case only one or more lens models are specified, the predicted image 
       positions are estimated by calculating the derivatives of :math:`\vec{\beta}(\vec{\theta})`
       and using them to compensate for small differences in source plane
       positions :math:`\vec{\beta}`; or if `useFSolve` is ``True``, several
       iterations will be used to approximate the true solution using scipy's
       `fsolve <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html>`_
       routine.

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

     - `useAverageBeta`: This can be either a boolean or a list. If ``True``, then the back-projected image points are averaged
       to estimate a single source position. If ``False``, then each back-projected image
       is used as an estimate of the source position. If this is a list, the number of
       entries should match `imgList`, and each entry is used as the source position for
       that `imgList` entry.

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

    if type(useAverageBeta) != bool:
        if len(useAverageBeta) != len(imgList):
            raise Exception("For predefined source plane positions, the number of source position should match the imgList length")

    cosmology, lensPlane, origLensModel, useTrace, createImgPlaneFn, allPoints = _commonInitFindOptRetrace(imgList, lensModel, cosmology, reduceImages)

    sources = []
    for srcIdx, (imgPos, z) in enumerate(allPoints):

        sourceInfo = [ ]
        if type(useAverageBeta) != bool:
            betas = [ useAverageBeta[srcIdx] ]
        else:
            betas = [ np.average([d["beta"] for d in imgPos], 0) ] if useAverageBeta else [ d["beta"] for d in imgPos ]

        if useTrace:
            imgPlane = createImgPlaneFn(lensPlane, z)

            observedThetas = [ x["theta"] for x in imgPos ]
            predictedThetas = [ [] for i in range(len(observedThetas)) ]
            betaFromThetaPred = [ [] for i in range(len(observedThetas)) ]
            # Trace each beta, and align the resulting position with the observations
            for beta in betas:
                thetas = imgPlane.traceBetaApproximately(beta)
                thetas = _align(thetas, observedThetas, maxPermSize)

                for i in range(len(thetas)):
                    predictedThetas[i].append(thetas[i])
                    betaFromThetaPred[i].append(imgPlane.traceTheta(thetas[i]))

            for i in range(len(observedThetas)):
                sourceInfo.append({ "theta_obs": observedThetas[i],
                                    "theta_pred": predictedThetas[i],
                                    "beta_from_theta_pred": betaFromThetaPred[i],
                                    "beta_est": betas})
        else:
            for x in imgPos:
                if useFSolve:
                    imgPlane = createImgPlaneFn(lensPlane, x["z"])

                thetaPred = [ ]
                betaFromThetaPred = []
                theta0, beta0, invderiv0 = x["theta"], x["beta"], x["invbetaderivs"]
                for beta1 in betas:
                    if useFSolve:
                        # We're looking for the images of beta1, starting from the
                        # estimate at observed theta0
                        th = _findRealTheta(imgPlane, beta1, theta0)
                    else:
                        dbeta = beta1-beta0
                        dtheta = invderiv0.dot(dbeta)
                        th = theta0 + dtheta

                    thetaPred.append(th)
                    betaFromThetaPred.append(imgPlane.traceTheta(th))

                sourceInfo.append({ "theta_obs": theta0,
                                    "theta_pred": thetaPred,
                                    "beta_from_theta_pred": betaFromThetaPred,
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


def createMonopoleBasisFunctions(avoidSources, Dd, subDiv, size, center = [0, 0], 
                                  widthFactor = 3.0, rangeFactor = 4.0,
                                  centralDensity = 1.0,
                                  overlapNeeded = None,
                                  randomOffset = True,
                                  cellCenterCallback = None,
                                  cellCenterCallbackState = None):
    """This creates a set of the monopole basis functions from
    `2008MNRAS.389..415L <https://ui.adsabs.harvard.edu/abs/2008MNRAS.389..415L/abstract>`_
    with which mass can be redistributed in between the images in a lensing system.
    The use from this article is illustrated in the `example_monopoledegen <https://github.com/j0r1/GRALE2/tree/master/inversion_examples/example_monopoledegen>`_
    example.

    It is also used in a different way in `2012MNRAS.425.1772L <https://ui.adsabs.harvard.edu/abs/2012MNRAS.425.1772L/abstract>`_, 
    where the deflection field at some image points is modified. This kind of use is
    demonstrated in example `example_gen_msd <https://github.com/j0r1/GRALE2/tree/master/inversion_examples/example_gen_msd>`_.

    Arguments:
     - `avoidSources`: a list of :class:`ImagesData <grale.images.ImagesData>` instances,
       describing the images to avoid. None of the generated monopole basis functions will
       overlap with these images.
     - `Dd`: the angular diameter distance to the lens plane, to be used in these basis 
       functions.
     - `subDiv`, `size`, `center`: the square region of size `size`, centered on `center`,
       will be subdivided into `subDiv` cells along each axis.
     - `widthFactor`: if a monopole can be added for a certain grid cell, the width of
       the positive part will be this factor times the cell width.
     - `rangeFactor`: for each cell, the minimal distance to the images will be determined.
       Only if this is larger than this factor times the cell width, will a monopole be
       added.
     - `centralDensity`: the monopole basis function will have this central density.
     - `overlapNeeded`: in the case that we're trying to modify properties at some points,
       we need to make sure that a basis function actually covers these points, otherwise
       it will not have any effect. This parameter is a set of :class:`ImagesData <grale.images.ImagesData>`
       instances, with which overlap is required.
     - `randomOffset`: if set to ``True`` (the default), the specified center will not be
       used exactly, but some randomness will be added.
     - `cellCenterCallback`: if set, for each cell of the resulting grid this callback
       function will be called. It will be passed four parameters: the center of the cell,
       the size of the cell, a boolean indicating if a monopole basis function could be
       added for this position, and the `cellCenterCallbackState` argument that's specified
       next.
     - `cellCenterCallbackState`: when the previous callback function is set, this value
       will be used as it's last argument.
    """
    
    import grale.lenses as lenses
    import grale.grid as ggrid
    import random

    if overlapNeeded:
        overlapNeeded = np.array([ pt["position"] for img in overlapNeeded for i in img.getAllImagePoints() for pt in i])
    else:
        overlapNeeded = None
    
    def getBorder(s, i):
        try:
            return s.getBorder(i)
        except:
            pass
        
        return s.getConvexHull(i)

    borders = [ [ getBorder(s, i) for i in range(s.getNumberOfImages()) ] for s in avoidSources ]

    def getDistanceToImages(pt):
        from shapely.geometry import Point, Polygon

        minDistSquared = float("inf")
        minPt = None
        for srcBorders in borders:
            for border in srcBorders:
                if Point(pt).within(Polygon(border)): # Check if point is inside polygon
                    return 0.0, pt.copy()

                for imgPt in border:
                    diff = imgPt - pt
                    distSquared = diff[0]**2 + diff[1]**2
                    if distSquared < minDistSquared:
                        minPt = imgPt
                        minDistSquared = distSquared
                        
        return minDistSquared**0.5, minPt

    cellSize = size/subDiv
    if randomOffset:
        rndX = (random.random()-0.5)*cellSize
        rndY = (random.random()-0.5)*cellSize
    else:
        rndX, rndY = 0, 0
    
    g = ggrid.createUniformGrid(size, np.array(center, dtype=np.double) + np.array([rndX, rndY]), subDiv)
    g = ggrid.fractionalGridToRealGrid(g)
    
    basisFunctions = []
    
    factor = widthFactor
    
    for cell in g:
        ctr = np.array(cell["center"])
        hasBasisFunction = False

        minDist, minPt = getDistanceToImages(ctr)
        minDist *= 0.99
        
        if minDist > rangeFactor*cellSize:
            
            if overlapNeeded is None:
                found = True
            else:
                # check if it overlaps with an image point, if not, there's no point in adding
                # the monopole since it won't affect the deflection angle at any of the image
                # points
                diffs = overlapNeeded - ctr
                lengths = np.sum(diffs*diffs, 1)**0.5
                isSmaller = lengths < minDist
                found = np.any(isSmaller)                
                
            if found:
                zeroPoint = (cellSize*factor)/minDist
                mp = lenses.ZeroMassLens(Dd, { 
                    "density": centralDensity,
                    "radius": minDist,
                    "zeropoint": zeroPoint
                })

                basisFunctions.append({
                    "lens": mp,
                    "center": ctr,
                    "mass": 1 # This is not used, but cannot be set to zero
                })
                hasBasisFunction = True

        if cellCenterCallback:
            cellCenterCallback(ctr, cellSize, hasBasisFunction, cellCenterCallbackState)
        
    return basisFunctions

def createThetaGrid(bottomLeft, topRight, NX, NY):
    """This creates a grid of xy-position, from `bottomLeft` to `topRight`,
    where the x-direction is evenly divided into `NX` points, and the
    y-direction into `NY` points. The result is a numpy array of
    shape (NY, NX, 2)"""

    thetas = np.empty([NY,NX,2], dtype=np.double)
    thetas[:,:,0], thetas[:,:,1] = np.meshgrid(np.linspace(bottomLeft[0], topRight[0], NX), 
                                               np.linspace(bottomLeft[1], topRight[1], NY))
    return thetas

def createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, maskRegions, enlargements=2, 
                                 enlargeDiagonally=False, circleToPolygonPoints=10000):
    """The `bottomLeft`, `topRight`, `NX` and `NY` parameters are passed to
    :func:`createThetaGrid <grale.util.createThetaGrid>`. Then, for this grid a mask
    is created based on the other arguments.

    The `maskRegions` argument is a list of region descriptions that are combined
    into one final mask for the grid, a boolean numpy array of shape (NY, NX). The
    entries of the list can be the following:

     - an :class:`ImagesData <grale.images.ImagesData>` instance: the grid positions
       covered by the images in this instance are set to ``True`` in the mask.
     - a 2D boolean array, also of shape (NY, NX): the positions where this mask
       is ``True`` are also set in the final mask.
     - a dictionary, with at least a key called ``type``, which can take on the
       following values:

        - ``hull``: in this case, a key ``imgdata`` should be present as well. It should
          contain either a single :class:`ImagesData <grale.images.ImagesData>` instance,
          or a list of these instances. All points will be collected, and the convex
          hull of this region will be used.
        - ``circle``: for this type, ``center`` and ``radius`` should be present as well.
          Internally this will be converted into a polygon, of which the number of
          points can be controlled by the `circleToPolygonPoints` argument.
        - ``rect``: describes a rectangle, and ``bottomleft`` and ``topright`` properties
          should also be present.
        - ``square``: for this type, ``center`` and ``size`` should also be present.
        - ``polygon``: this describes a polygon, there should be a ``coord`` property
          that represents a list of points that make up the polygon.
        - ``point``: for this type, a ``coord`` property describes a single point.

        For each of these types, an ``invert`` key may be present as well. If set to
        ``True``, then the region will be inverted.

    The `enlargements` value determines how many times the mask should grow. If a
    value is ``True``, then the values to the left, right, below and above will
    be set to ``True`` as well. If `enlargeDiagonally` is set to ``True``, the
    four diagonal values will also be set to ``True``. This procedure is repeated
    the number of times specified by `enlargements`.

    The return value is a tuple of the resulting grid, and the final mask.
    """
    from . import images
    import numpy as np

    V = lambda x, y: np.array([x, y], dtype=np.double)

    def extractImageRegion(img, imgIdx, tryBorder = True):
        numPoints = img.getNumberOfImagePoints(imgIdx)
        if numPoints == 0:
            print("WARNING: ignoring image idx {} with no points".format(imgIdx))
            return []

        if numPoints <= 3:
            return [ img.getImagePointPosition(imgIdx, pt) for pt in range(numPoints) ]

        # Try border
        if tryBorder and img.hasTriangulation():
            try:
                return img.getBorder(imgIdx)
            except images.ImagesDataException: # Assuming no border found
                return img.getConvexHull(imgIdx)

        # Use convex hull
        return img.getConvexHull(imgIdx)

    def extractImageRegions(img):
        l = []
        for i in range(img.getNumberOfImages()):
            l.append({ "type": "polygon", "coord": extractImageRegion(img, i), "invert": False })
        return l

    newRegionList = []
    updatedRegionList = False
    for r in maskRegions:
        if type(r) == images.ImagesData:
            updatedRegionList = True
            newRegionList += extractImageRegions(r)

        elif isinstance(r, np.ndarray):
            newRegionList.append(r)

        elif r["type"] == "hull":
            
            updatedRegionList = True
            
            allPts = images.ImagesData(1)
            
            l = [ r["imgdata"] ] if type(r["imgdata"]) is not list else r["imgdata"]
            for imgDat in l:
                for img in range(imgDat.getNumberOfImages()):
                    for pt in range(imgDat.getNumberOfImagePoints(img)):
                        xy = imgDat.getImagePointPosition(img, pt)
                        allPts.addPoint(0, xy)
            
            newRegionList.append({
                "type": "polygon",
                "coord": allPts.getConvexHull(0),
                "invert": r["invert"] if "invert" in r else False
            })
        elif r["type"] == "circle": # Convert to polygon
            
            updatedRegionList = True
            
            ctr = r["center"]
            radius = r["radius"]
            
            numPoints = circleToPolygonPoints + 1
            
            newRegionList.append({
                "type": "polygon",
                "coord": [ ctr + radius*V(np.cos(angle), np.sin(angle)) for angle in np.linspace(0, 2*np.pi, numPoints) ],
                "invert": r["invert"] if "invert" in r else False
            })
        elif r["type"] == "rect": # Convert to polygon
            updatedRegionList = True
            bl, tr = r["bottomleft"], r["topright"]
            x0, y0 = bl
            x1, y1 = tr
            newRegionList.append({
                "type": "polygon",
                "coord": [ V(x0, y0), V(x1, y0), V(x1, y1), V(x0, y1), V(x0, y0) ],
                "invert": r["invert"] if "invert" in r else False
            })
        elif r["type"] == "square": # Convert to polygon
            updatedRegionList = True
            ctr = r["center"]
            sz = r["size"]
            bl = ctr-sz/2
            tr = ctr+sz/2
            x0, y0 = bl
            x1, y1 = tr
            newRegionList.append({
                "type": "polygon",
                "coord": [ V(x0, y0), V(x1, y0), V(x1, y1), V(x0, y1), V(x0, y0) ],
                "invert": r["invert"] if "invert" in r else False
            })

        else:
            if not "invert" in r:
                r2 = r.copy()
                r2["invert"] = False
            else:
                r2 = r

            newRegionList.append(r2)
            
    if updatedRegionList:
        return createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, newRegionList, enlargements, enlargeDiagonally,
                                            circleToPolygonPoints)
    
    thetas = createThetaGrid(bottomLeft, topRight, NX, NY)
    
    import cairo
    
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, NX, NY)
    ctx = cairo.Context(surface)
    ctx.set_antialias(cairo.ANTIALIAS_NONE)
    ctx.scale(1, 1)
    
    ctx.set_source_rgb(0, 0, 0)
    ctx.paint()
    
    def translateCoord(xy):
        return (xy - bottomLeft)/(topRight-bottomLeft) * V(NX-1, NY-1) + V(0.5, 0.5)

    ctx.set_source_rgb(1, 1, 1)

    mask = np.zeros((NY,NX), dtype=bool)
    
    for region in maskRegions:

        # Check if it's a numpy array itself
        if isinstance(region, np.ndarray):
            if len(region.shape) != 2 or region.shape[0] != NY or region.shape[1] != NX:
                raise Exception("Region specified as ndarray is not shape compatible with the mask: expecting ({},{}), got {}".format(NY, NX, region.shape))

            region = region.astype(bool)
            mask |= region
            continue

        # Not an array, need to draw something and interpret that
        tp = region["type"]
        if tp == "polygon":
            points = [ translateCoord(c) for c in region["coord"] ]
            ctx.set_line_width(1)

            ctx.move_to(*points[0])
            for xy in points[1:]:
                ctx.line_to(*xy)
            ctx.close_path()
            ctx.fill()
            
            ctx.move_to(*points[0])
            for xy in points[1:]:
                ctx.line_to(*xy)
            ctx.close_path()
            ctx.stroke()
            
            for xy in points:
                ctx.rectangle(xy[0]-0.5, xy[1]-0.5, 1, 1)
                ctx.fill()
            
        elif tp == "point":

            xy = translateCoord(region["coord"])
            ctx.rectangle(xy[0]-0.5, xy[1]-0.5, 1, 1)
            ctx.fill()

        else:
            raise Exception("Unknown region type '{}''".format(tp))

        pixel_data = np.frombuffer(surface.get_data(), np.uint8)
        pixel_data = pixel_data.reshape((NY, NX, 4))
        regMask = pixel_data[:,:,0] > 0

        if region["invert"]:
            regMask = ~regMask

        mask |= regMask

    def shiftMask(m, dx, dy):
        NY, NX = m.shape
        m2 = np.zeros(m.shape, dtype=m.dtype)

        leftx0, leftx1, rightx0, rightx1 = (dx,NX,0,NX-dx) if dx > 0 else (0,NX+dx,-dx,NX)
        lefty0, lefty1, righty0, righty1 = (dy,NY,0,NY-dy) if dy > 0 else (0,NY+dy,-dy,NY)
        m2[lefty0:lefty1,leftx0:leftx1] = m[righty0:righty1,rightx0:rightx1]
        return m2

    def enlargeMask(m, enlargeOffsets):
        mNew = m.copy()
        for dx,dy in enlargeOffsets:
            mNew |= shiftMask(m, dx, dy)
        return mNew

    enlargeOffsets = [ (-1,0), (1,0),(0,-1),(0,1), ]
    if enlargeDiagonally:
        enlargeOffsets += [ (-1,-1),(-1,1),(1,-1),(1,1) ]
    
    for i in range(enlargements):
        mask = enlargeMask(mask, enlargeOffsets)
    
    return thetas, mask
