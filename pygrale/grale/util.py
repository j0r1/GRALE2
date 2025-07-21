"""Module meant for various utilities."""

import numpy as np
from numpy.linalg import inv
from .inverters import _getNumHelpers
from .images import ImagesData
from .constants import ANGLE_ARCSEC, MASS_SUN
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

_useOldFindRealTheta = False

def localretrace_fsolve(*args):

    if len(args) == 2:
        return None

    imgPlane, beta0, theta, unusedOptions = args

    from .constants import ANGLE_ARCSEC
    import scipy.optimize as opt

    #startDiff = (imgPlane.traceTheta(theta) - beta)/ANGLE_ARCSEC

    cachedValues = [ None, None, None ]
    def getBetaAndDerivs(x):
        if cachedValues[0] is not None and x[0] == cachedValues[0][0] and x[1] == cachedValues[0][1]:
            beta, derivs = cachedValues[1:]
        else:
            beta, derivs = imgPlane.getBetaAndDerivatives(x)
            cachedValues[0] = x.copy()
            cachedValues[1] = beta.copy()
            cachedValues[2] = derivs.copy()
        return beta, derivs

    def f(x):
        beta, derivs = getBetaAndDerivs(x)
        return (beta - beta0)/ANGLE_ARCSEC

    def fprime(x):
        beta, derivs = getBetaAndDerivs(x)
        return derivs

    def f_old(x):
        return (imgPlane.traceTheta(x) - beta0)/ANGLE_ARCSEC

    if _useOldFindRealTheta:
        r = opt.fsolve(f, x0=theta)
    else:
        r = opt.fsolve(f, fprime=fprime, x0=theta)

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
            beta = imgPlane.traceTheta(theta)
            imgPos.append({ "theta": theta, "beta": beta, "z": z })

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

def _processBetas_noreduction(imgPlane, srcIdx, imgPos):
    return [ d["beta"] for d in imgPos ]

def _processBetas_average_unweighted(imgPlane, srcIdx, imgPos):
    return [ np.average([d["beta"] for d in imgPos], 0) ]

def _processBetas_average_weighted_mu(imgPlane, srcIdx, imgPos):
    totMag = 0
    avgBeta = np.array([0.0, 0.0])
    
    for x in imgPos:
        beta, derivs = imgPlane.getBetaAndDerivatives(x["theta"])
        assert np.sum((beta - x["beta"])**2) == 0 # Sanity check

        mag = 1.0/abs(derivs[0,0]*derivs[1,1] - derivs[1,0]*derivs[0,1])
        totMag += mag
        avgBeta += beta*mag

    avgBeta /= totMag
    return [ avgBeta ]

def calculateImagePredictions(imgList, lensModel, cosmology=None,
                     reduceImages="average",
                     reduceSources="average", # or "magweighted", or "noreduction", or list, or function
                     maxPermSize=7,
                     localTraceFunction = "fsolve",
                     localTraceFunctionOptions = None,
                     useFSolve = None, # deprecated
                     useAverageBeta=None, # deprecated
                    ):
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

     - `reduceImages`: each image in an images data set will be converted to a single
       point. By default, this is done by averaging the image point positions. If something
       else needs to be done, you can specify a function here, which will be called with
       two arguments: the :class:`ImagesData <grale.images.ImagesData>` instance and the
       image index corresponding to the particular image.

     - `reduceSources`: TODO "average", "magweighted", "noreduction", list, or function

     - `maxPermSize`: when a source plane position is used to estimate the corresponding
       image plane positions, these predicted positions are grouped with the observed
       positions. If the derivatives-based method is used, then this is automatically
       possible, but if the more accurate tracing is used, we need to find out which
       predicted point corresponds to which observed point. As a first attempt, the
       procedure tries to use the points closest to the observed ones, but this may not
       be possible and different permutations of the positions will then be explored.
       Since this can become quite slow, an error is generated if more elements than
       this would need to be permutated.

     - `localTraceFunction`: TODO

     - `localTraceFunctionOptions`: TODO

    """
    # - `useFSolve`: deprecated now
    # if only a lens model is used, this cause SciPy's `fsolve <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.optimize.fsolve.html>`_
    #   to be used to zoom in on the real image plane vectors. Otherwise, a one-step
    #   approach is used based on the deflection angle derivatives.
    #- `useAverageBeta`: TODO: deprecated
    #   This can be either a boolean or a list. If ``True``, then the back-projected image points are averaged
    #   to estimate a single source position. If ``False``, then each back-projected image
    #   is used as an estimate of the source position. If this is a list, the number of
    #   entries should match `imgList`, and each entry is used as the source position for
    #   that `imgList` entry.

    # localTraceFunction first called with allPoints, localTraceFunctionOptions
    # then with (lensObj, traceTheta, getBetaAndDerivatives, returnedOptions)
    # and should return optimized theta

    if useFSolve is not None:
        raise Exception("'useFSolve' parameter is deprecated, use 'localTraceFunction'")
    if useAverageBeta is not None:
        raise Exception("'useAverageBeta' parameter is deprecated, use 'reduceSources'")

    # TODO: use each beta as a possible source
    reduceSourcesFunction = None
    if type(reduceSources) == str:
        if reduceSources == "noreduction":
            reduceSourcesFunction = _processBetas_noreduction
        elif reduceSources == "average":
            reduceSourcesFunction = _processBetas_average_unweighted
        elif reduceSources == "magweighted":
            reduceSourcesFunction = _processBetas_average_weighted_mu
        else:
            raise Exception(f"Unknown value '{reduceSources}' for 'reduceSources' parameter, expecting 'noreduction', 'average' or 'magweighted'")

    elif type(reduceSources) == list:
        reduceSourcesFunction = lambda imgPlane,srcIdx,imgPos: reduceSources[srcIdx] # Just return part of the list

    else:
        reduceSourcesFunction = reduceSources # This parameter specifies the function to use

    #if type(useAverageBeta) != bool:
    #    if len(useAverageBeta) != len(imgList):
    #        raise Exception("For predefined source plane positions, the number of source position should match the imgList length")

    cosmology, lensPlane, origLensModel, useTrace, createImgPlaneFn, allPoints = _commonInitFindOptRetrace(imgList, lensModel, cosmology, reduceImages)

    if not useTrace: # Not using a rasterized image plane, use local function

        if localTraceFunction == "fsolve":
            if localTraceFunctionOptions is not None:
                raise Exception("'localTraceFunctionOptions' should be None for 'fsolve'")
            localTraceFunction = localretrace_fsolve

        elif type(localTraceFunction) == str and (localTraceFunction == "pgraleretrace" or localTraceFunction.startswith("pgraleretrace,")):
            # First get the correct functions, then call them
            localTraceFunction, localTraceFunctionOptions = _processPGraleRetraceType(localTraceFunction, localTraceFunctionOptions)
            localTraceFunctionOptions = localTraceFunction(allPoints, localTraceFunctionOptions)

        # TODO: other default solvers?

        else: # Assume it's a function
            localTraceFunctionOptions = localTraceFunction(allPoints, localTraceFunctionOptions)

    sources = []
    for srcIdx, (imgPos, z) in enumerate(allPoints):

        imgPlane = createImgPlaneFn(lensPlane, z)

        sourceInfo = [ ]

        betas = reduceSourcesFunction(imgPlane, srcIdx, imgPos)
        assert len(betas) == 1 or len(betas) == len(imgPos)

        if useTrace:

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

                thetaPred = [ ]
                betaFromThetaPred = []
                theta0, beta0 = x["theta"], x["beta"]
                for beta1 in betas:
                    # We're looking for the images of beta1, starting from the
                    # estimate at observed theta0
                    th = localTraceFunction(imgPlane, beta1, theta0, localTraceFunctionOptions)
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

def adjustShearMeasurements(pixelFrameCoords, pixelFrameGamma1, pixelFrameGamma2, centeredRaDecCoords, mirror=False, tol=None):
    """It is possible that shear measurements are specified in one coordinate
    frame, but you need to know them in a rotated, possibly mirrored frame. This
    function can perform the required calculation. In the names of the parameters
    we're assuming that shear measurements need to be transformed from values that
    are given in a coordinate frame that's aligned with the pixels of the CCD camera,
    and need to be transformed to a coordinate system that's aligned with RA/Dec
    coordinates (but also recentered). Of course, these names just suggest coordinate
    systems, the transformation could be the other way around as well, for example.

    Arguments:
     - `pixelFrameCoords`: The coordinates of the shear measurements in the coordinate
       frame that's aligned with the CCD pixels.
     - `pixelFrameGamma1`: The first shear component as measured in the coordinate
       frame that's aligned with the pixels.
     - `pixelFrameGamma2`: Same for the second shear component.
     - `centeredRaDecCoords`: The same points as specified in `pixelFrameCoords`, but
       now having their coordinates in the frame that's being transformed to, for example
       one aligned with the RA/Dec axes.
     - `mirror`: In case the two coordinate systems don't differ by simply a rotation,
       this flag can be set to indicate this.
     - `tol`: tolerance parameter that is passed to SciPy's minimize function ("Nelder-Mead")
       is used.
    """
       
    origShape1 = pixelFrameGamma1.shape
    origShape2 = pixelFrameGamma2.shape
    pixelFrameCoords = pixelFrameCoords.copy().reshape((-1,2))
    centeredRaDecCoords = centeredRaDecCoords.copy().reshape((-1,2))
    
    if mirror:
        centeredRaDecCoords[:,0] = -centeredRaDecCoords[:,0]
    
    def getScale(coords):
        return np.sum((np.max(coords,axis=0) - np.min(coords, axis=0))**2)**0.5
    
    pfScale = getScale(pixelFrameCoords)
    pixelFrameCoords /= pfScale
    rdScale = getScale(centeredRaDecCoords)
    centeredRaDecCoords /= rdScale
    
    startOffset = np.mean(centeredRaDecCoords, axis=0) - np.mean(pixelFrameCoords, axis=0)
    
    # Find rotation angle
    def getRotMatrix(rotAngle):
        return np.array([
            [ np.cos(rotAngle), -np.sin(rotAngle) ],
            [ np.sin(rotAngle), np.cos(rotAngle) ]
        ])
    
    def fitness(params):
        x0, y0, scale, rotAngle = params
        rotScaledOffset = np.matmul(getRotMatrix(rotAngle), pixelFrameCoords.T).T*scale
        rotScaledOffset[:,0] += x0
        rotScaledOffset[:,1] += y0
        
        return np.sum((centeredRaDecCoords - rotScaledOffset)**2)
        
    from scipy.optimize import minimize
    r = minimize(fitness, (startOffset[0], startOffset[1], 1, 0), method="Nelder-Mead", tol=tol)
    if not r.success:
        raise Exception("Couldn't optimize to find rotation angle between frames")
    
    angle = r.x[3]
    
    pixelFrameGamma = np.concatenate((pixelFrameGamma1.reshape((-1,1)), pixelFrameGamma2.reshape(-1,1)), axis=1)
    rotGamma = np.matmul(getRotMatrix(2*angle), pixelFrameGamma.T).T
    
    if mirror:
        rotGamma[:,1] = -rotGamma[:,1]
    return rotGamma[:,0].reshape(origShape1), rotGamma[:,1].reshape(origShape2), angle

def pgraleretrace_noTrace(*args):
    if len(args) == 2:
        return args[1]

    ipObj, betaTgt, thetaStart, opts = args

    if opts:
        _checkKeys(opts, ["type"], "Unexpected option for 'NoTrace'")

    return thetaStart

def pgraleretrace_singleStepNewton(*args):
    if len(args) == 2:
        return args[1]

    ipObj, betaTgt, thetaStart, opts = args

    if opts:
        _checkKeys(opts, ["type"], "Unexpected option for 'SingleStepNewton'")

    beta, derivs = ipObj.getBetaAndDerivatives(thetaStart)
    invderivs = inv(derivs)
    dbeta = betaTgt - beta
    dtheta = invderivs@dbeta
    return thetaStart + dtheta

def pgraleretrace_multiStepNewton(*args):

    if len(args) == 2:
        _, opts = args
        if not opts:
            opts = { }

        from .inversion import getDefaultsForRetraceType
        defaults = getDefaultsForRetraceType("MultiStepNewton")
        for k in defaults:
            if not k in opts:
                opts[k] = defaults[k]
        return opts # object containing number of steps

    ipObj, betaTgt, thetaStart, opts = args

    _checkKeys(opts, [ "type", "evaluations"],
               "Unexpected key for 'MultiStepNewton' options")

    numSteps = opts["evaluations"]

    return _multiStepNewton_singlerun(ipObj, thetaStart, betaTgt, numSteps)[0]

def _multiStepNewton_singlerun(ipObj, thetaStart, betaTarget, numIterations):

    theta = thetaStart
    betaCur, betaDerivs = ipObj.getBetaAndDerivatives(theta)
    betaDiff = betaTarget - betaCur

    bestBetaDiffSize = (betaDiff[0]**2 + betaDiff[1]**2)**0.5
    bestRetraceTheta = theta

    while numIterations > 0:
        thetaDiff = inv(betaDerivs)@betaDiff
        #print("thetaDiff =", thetaDiff/ANGLE_ARCSEC)

        while numIterations > 0:
            numIterations -= 1

            betaCur, betaDerivs = ipObj.getBetaAndDerivatives(theta + thetaDiff)
            betaDiff = betaTarget - betaCur

            betaDiffSize = (betaDiff[0]**2 + betaDiff[1]**2)**0.5
            if betaDiffSize < bestBetaDiffSize:
                theta = theta + thetaDiff
                bestBetaDiffSize = betaDiffSize
                bestRetraceTheta = theta
                break

            else:
                thetaDiff *= 0.5

    return bestRetraceTheta, bestBetaDiffSize

def _findRetraceTheta_level(ipObj, levelCoords, dxy, thetaStart, betaTarget, acceptThreshold, numIterations):

    totalBestBetaDiff = np.inf
    totalBestRetraceTheta = np.array([np.inf, np.inf], dtype=np.float64)

    closestAcceptableThetaDiff2 = np.inf
    closestAcceptableBetaDiff = np.inf
    closestBestRetraceTheta = np.array([np.inf, np.inf], dtype=np.float64)

    for DX, DY in levelCoords:
        dt = np.array([DX, DY], dtype=np.float64)*dxy
        theta = thetaStart + dt

        curBestRetraceTheta, curBestBetaDiff = _multiStepNewton_singlerun(ipObj, theta, betaTarget, numIterations)
        if curBestBetaDiff < totalBestBetaDiff:
            totalBestBetaDiff = curBestBetaDiff
            totalBestRetraceTheta = curBestRetraceTheta

        if curBestBetaDiff <= acceptThreshold:
            retrDiff = curBestRetraceTheta - thetaStart
            thetaDiff2 = retrDiff[0]**2 + retrDiff[1]**2
            if thetaDiff2 < closestAcceptableThetaDiff2:
                closestAcceptableBetaDiff = curBestBetaDiff
                closestAcceptableThetaDiff2 = thetaDiff2
                closestBestRetraceTheta = curBestRetraceTheta

    if closestAcceptableThetaDiff2 < np.inf:
        return closestBestRetraceTheta, closestAcceptableBetaDiff

    return totalBestRetraceTheta, totalBestBetaDiff

def pgraleretrace_expandedMultiStepNewton(*args):
    if len(args) == 2:
        allImages, opts = args

        from .inversion import getDefaultsForRetraceType
        defaultOpts = getDefaultsForRetraceType("ExpandedMultiStepNewton")
        if opts is None:
            opts = defaultOpts
        else:
            # Merge options
            opts = opts.copy()
            for k in defaultOpts:
                if not k in opts:
                    opts[k] = defaultOpts[k]

        imgList = []
        for imgPoints, z in allImages:
            imgDat = ImagesData(len(imgPoints))
            for idx, x in enumerate(imgPoints):
                imgDat.addPoint(idx, x["theta"])

            imgList.append({
                "images": imgDat,
                "params": { "type": "bayesstronglensing" } # Needed to use another internal function
            })

        from .inversion import _getAcceptThresholdFromBayesStrongLensingImages, _getGridSpacingFromBayesStrongLensingImages
        if opts["acceptthreshold"] == "auto":
            opts["acceptthreshold"] = _getAcceptThresholdFromBayesStrongLensingImages(imgList)

        if opts["gridspacing"] == "auto":
            opts["gridspacing"] = _getGridSpacingFromBayesStrongLensingImages(imgList)

        return opts

    ipObj, betaTgt, thetaStart, opts = args

    _checkKeys(opts, ["type", "maxgridsteps", "evalsperpoint", "acceptthreshold", "gridspacing", "layout"],
               "Unexpected key in ExpandedMultiStepNewton options")

    from .inversionparams import getCoordinatesForExpandedMultiStepNewtonGridSteps
    stepInfo = getCoordinatesForExpandedMultiStepNewtonGridSteps(opts)

    assert opts["maxgridsteps"] == len(stepInfo)
    evalsPerPoint = opts["evalsperpoint"]
    acceptThreshold = opts["acceptthreshold"]
    dxy = opts["gridspacing"]

    totalBestBetaDiff = np.inf
    totalBestRetraceTheta = np.array([np.inf, np.inf], dtype=np.float64)

    for idx, levelCoords in enumerate(stepInfo):
        #print("processing level", idx)
        curBestRetraceTheta, curBestBetaDiffSize = _findRetraceTheta_level(ipObj, levelCoords, dxy, thetaStart,
                                                                          betaTgt, acceptThreshold, evalsPerPoint)

        if curBestBetaDiffSize <= acceptThreshold:
            #print("Found at level", idx, curBestBetaDiffSize/ANGLE_ARCSEC)
            return curBestRetraceTheta

        if curBestBetaDiffSize < totalBestBetaDiff:
            totalBestBetaDiff = curBestBetaDiffSize
            totalBestRetraceTheta = curBestRetraceTheta

    #print("Not found, best is", totalBestBetaDiff/ANGLE_ARCSEC)
    return totalBestRetraceTheta

def _processPGraleRetraceType(name, localTraceFunctionOptions):
    nameParts = name.split(",")
    name = nameParts[0]
    if name != "pgraleretrace":
        raise Exception("Unexpected: algorithl name should be 'pgraleretrace'")

    if len(nameParts) == 1:
        opts = { "type": "ExpandedMultiStepNewton" } # Same default as in inversion
    else:
        algType = ",".join(nameParts[1:])
        opts = { "type": algType }
        if localTraceFunctionOptions and "type" in localTraceFunctionOptions:
            if localTraceFunctionOptions["type"] != algType:
                raise Exception(f"Type '{algType}' from short algorithm name does not match the one in 'localTraceFunctionOptions'")

    if localTraceFunctionOptions:
        for k in localTraceFunctionOptions:
            opts[k] = localTraceFunctionOptions[k]

    optsType = opts["type"]
    if optsType == "NoTrace":
        alg = pgraleretrace_noTrace
    elif optsType == "SingleStepNewton":
        alg = pgraleretrace_singleStepNewton
    elif optsType == "MultiStepNewton":
        alg = pgraleretrace_multiStepNewton
    elif optsType == "ExpandedMultiStepNewton":
        alg = pgraleretrace_expandedMultiStepNewton
    else:
        raise Exception(f"Unknown algorithm name '{optsType}' for 'pgraleretrace'")

    return alg, opts

def _checkKeys(obj, keyNames, errMsg):
    for k in obj:
        if not k in keyNames:
            raise Exception(errMsg + ": " + k)

def readLenstoolBayesDat(fn, columnsToRemove = [ "Nsample", "ln(Lhood)", "Chi2" ]):
    """TODO"""
    colNames = []
    for l in open(fn, "rt"):
        if not l.startswith("#"):
            break
        colNames.append(l[1:].strip())

    import pandas as pd
    df = pd.read_csv(fn, comment='#', sep=r'\s+', header=None, names = colNames)
    for n in columnsToRemove:
        del df[n]
    return df

_defaultRenameParams = {
    "angle": (1.0, "degree"),
    "angdist": (ANGLE_ARCSEC, "arcsec"),
    "mass": (MASS_SUN, "msun"),
    "density": (1.0, "kgm2"),
    "veldisp": (1000.0, "kms"),
}

def _defaultRenameColumnsFunction(colName, params):
    origColName = colName
    
    lensNr = None
    lensPrefix = "lens_"
    anglePrefix = "angle_"
    factorPrefix = "factor_"
    
    def _processType(colName, typeName):
        scale, suffix = params[typeName]
        newColName = colName if suffix is None else colName + "_" + suffix
        return newColName, scale
        
    if colName.startswith(anglePrefix):
        nr = int(colName[len(anglePrefix):]) # Just a check to see if it parses, not really needed here
        return _processType(colName, "angle")
    
    if colName.startswith(factorPrefix):
        nr = int(colName[len(factorPrefix):]) # Again, just to see if it parses
        return colName, 1.0
    
    # This piece of code determines a lens number, and then updates the column name
    # to what's after the comma
    if colName.startswith(lensPrefix):
        idx = colName.find(",")
        assert idx > len(lensPrefix), f"Couldn't find ',' in column name that starts with 'lens_': {colName}"
        lensNr = int(colName[len(lensPrefix):idx])
        colName = colName[idx+1:]

    # Remove the '_scaled' suffix if present
    scaledSuffix = "_scaled"
    if colName.endswith(scaledSuffix):
        colName = colName[:-len(scaledSuffix)]
    
    # It's possible that this is just some property of a composite lens, or of a multiple
    # plummer lens
    if colName.startswith("x_") or colName.startswith("y_"):
        nr = int(colName[2:])
        return _processType(colName if lensNr is None else colName[:2] + f"{lensNr}-{nr}", "angdist")
    
    if colName.startswith("mass_"): # Component of multiple plummer lens
        nr = int(colName[5:])
        return _processType(colName if lensNr is None else colName[:5] + f"{lensNr}-{nr}", "mass")
    
    if colName.startswith("width_"): # Component of multiple plummer lens
        nr = int(colName[6:])
        return _processType(colName if lensNr is None else colName[:6] + f"{lensNr}-{nr}", "angdist")
    
    singleProperties = {
        "centraldensity": "density",
        "core": "angdist",
        "coreradius": "angdist",
        "density": "density",
        "ellipticity": None,
        "epsilon": None,
        "mass": "mass",
        "scaleradius": "angdist",
        "shearangle": None,
        "shearsize": None,
        "sigma": "veldisp",
        "velocitydispersion": "veldisp",
        "width": "angdist"
    }
    
    for n in singleProperties:
        if n == colName:
            colName = colName if lensNr is None else colName + f"_{lensNr}"
            unit = singleProperties[n]
            if unit is None:
                return colName, 1.0
            return _processType(colName, unit)
            
    raise Exception(f"Unexpected: unhandled column name {colName} (from {origColName})")

def readParametricSamples(parametersFileName, samplesFileName, renameColumnsFunction = _defaultRenameColumnsFunction, renameColumnsFunctionParams = _defaultRenameParams):
    """TODO"""

    # Read the information about the parameters that were being optimized. This was written
    # in the inversion script, and it contains a name for each parameter in the samples file
    # as well as a scale factor. This scale factor is needed to convert the values from the
    # samples file to regular units, as they use some arbitratry units to keep the calculated
    # values in a sensible range
    import pickle
    import pandas as pd

    paramInfo = pickle.load(open(parametersFileName, "rb"))
    names = [ x["name"] for x in paramInfo ]
    scales = np.array([ x["scalefactor"] for x in paramInfo ]).reshape((1,-1))

    # Read the samples file, reshape it so that it has some amount of rows, each containing the
    # values for the parameters in the sample. The amount of values per sample is ultimately
    # obtained from the "paraminfo.dat" file
    samples = np.fromfile(samplesFileName, dtype=np.float32).reshape((-1,len(names))).astype(np.float64)
    # Rescale the values to get the correct units
    scaledSamples = samples*scales

    df = pd.DataFrame(scaledSamples, columns=names).copy()
    if renameColumnsFunction is None:
        return df

    for n in names:
        r = renameColumnsFunction(n, renameColumnsFunctionParams)
        if r is None:
            continue

        oldName = n
        newName, scale = r
        df[newName] = df[oldName]/scale
        del df[oldName]

    return df
