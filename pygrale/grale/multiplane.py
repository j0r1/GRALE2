"""With the classes in this module you can simulate a multi-lensplane situation."""
from . import images
from . import contourfinder
from . import gridfunction
from . import privutil
import copy
import numpy as np

class MultiLensPlaneException(Exception):
    """An exception that's generated when something goes wrong in the classes in this module."""
    pass

class _MultiPlaneCache(object):
    def __init__(self, lz, thetas, renderer, feedbackObject, cosmology, calcDerivs):

        self.thetas = copy.deepcopy(thetas)
        self.lz = lz
        self.cosmology = cosmology
        self.calcDerivs = calcDerivs

        txxShape = thetas.shape[:-1]
        self.txxShape = txxShape

        # Calculate mappings of deflection angles and derivatives
        T, Txx, Txy, Tyx, Tyy = [], [], [], [], []
        alphas = { }
        derivs = { }
        alphasAndDerivs = { }

        def getAlphasAndDerivs(j):
            if j in alphasAndDerivs:
                ad = alphasAndDerivs[j]
            else:
                if renderer is None:
                    a = lz[j][0].getAlphaVector(T[j])
                    d = lz[j][0].getAlphaVectorDerivatives(T[j]) if calcDerivs else None
                else: # TODO: bypass derivatives if not requested?
                    result = renderer.renderXYVector(lz[j][0].toBytes(), T[j])
                    result = np.frombuffer(result, "double").reshape([-1,5]) # for each point: ax, ay, axx, ayy, axy
                    a = result[:,0:2].reshape(txxShape + (2,))
                    d = result[:,2:5].reshape(txxShape + (3,))

                ad = (a, d)
                alphasAndDerivs[j] = ad

            return ad

        def getAlphas(j):
            if j in alphas:
                return alphas[j]
            a,d = getAlphasAndDerivs(j)
            alphas[j] = a
            return a
        
        def getAlphaDeriv(j):
            if j in derivs:
                return derivs[j]

            a, tmp = getAlphasAndDerivs(j)
            d = [ tmp[:,:,0], tmp[:,:,2], tmp[:,:,2], tmp[:,:,1] ] 
            derivs[j] = d
            return d
        
        N = len(lz)
        for i in range(0, N):
            #print("i =", i)
            feedbackObject.onStatus("Processing lens {} at z = {:.2f}".format(i, lz[i][1]))

            Ti = copy.deepcopy(thetas)
            Di = cosmology.getAngularDiameterDistance(lz[i][1])
            
            if calcDerivs:
                Tixx = np.ones(txxShape, dtype=np.double)
                Tixy = np.zeros(txxShape, dtype=np.double)
                Tiyx = np.zeros(txxShape, dtype=np.double)
                Tiyy = np.ones(txxShape, dtype=np.double)
            
            for j in range(0, i):
                #print("j =", j)
                Dji = cosmology.getAngularDiameterDistance(lz[j][1], lz[i][1])
                Ti -= Dji/Di * getAlphas(j)
                if calcDerivs:
                    Tixx -= Dji/Di * ( getAlphaDeriv(j)[0] * Txx[j] + getAlphaDeriv(j)[1]*Tyx[j])
                    Tixy -= Dji/Di * ( getAlphaDeriv(j)[0] * Txy[j] + getAlphaDeriv(j)[1]*Tyy[j])
                    Tiyx -= Dji/Di * ( getAlphaDeriv(j)[2] * Txx[j] + getAlphaDeriv(j)[3]*Tyx[j])
                    Tiyy -= Dji/Di * ( getAlphaDeriv(j)[2] * Txy[j] + getAlphaDeriv(j)[3]*Tyy[j])
        
            T.append(Ti)
            if calcDerivs:
                Txx.append(Tixx)
                Txy.append(Tixy)
                Tyx.append(Tiyx)
                Tyy.append(Tiyy)

        # Make sure these are cached
        getAlphas(N-1)
        if calcDerivs:
            getAlphaDeriv(N-1)

        self.alphas = alphas
        self.derivs = derivs
        self.T = T
        self.Txx = Txx
        self.Tyy = Tyy
        self.Txy = Txy
        self.Tyx = Tyx

    def addSourcePlane(self, sourceRedshift):

        # Keep only the part that's relevant to us
        lz = self.lz[:]
        while lz and sourceRedshift < lz[-1][1]:
            del lz[-1]

        alphas = self.alphas
        derivs = self.derivs
        T = self.T
        Txx = self.Txx
        Tyy = self.Tyy
        Txy = self.Txy
        Tyx = self.Tyx

        cosmology = self.cosmology
        Ti = copy.deepcopy(self.thetas)
        Di = cosmology.getAngularDiameterDistance(sourceRedshift)
        if self.calcDerivs:
            Tixx = np.ones(self.txxShape, dtype=np.double)
            Tixy = np.zeros(self.txxShape, dtype=np.double)
            Tiyx = np.zeros(self.txxShape, dtype=np.double)
            Tiyy = np.ones(self.txxShape, dtype=np.double)
        
        N = len(lz)
        for j in range(0, N):
            Dji = cosmology.getAngularDiameterDistance(lz[j][1], sourceRedshift)
            Ti -= Dji/Di * alphas[j]
            if self.calcDerivs:
                Tixx -= Dji/Di * ( derivs[j][0] * Txx[j] + derivs[j][1]*Tyx[j])
                Tixy -= Dji/Di * ( derivs[j][0] * Txy[j] + derivs[j][1]*Tyy[j])
                Tiyx -= Dji/Di * ( derivs[j][2] * Txx[j] + derivs[j][3]*Tyx[j])
                Tiyy -= Dji/Di * ( derivs[j][2] * Txy[j] + derivs[j][3]*Tyy[j])

        if self.calcDerivs:
            return Ti, Tixx, Tixy, Tiyx, Tiyy
        return Ti

class MultiLensPlane(object):

    @staticmethod
    def _createThetaGrid(bottomLeft, topRight, numX, numY):
        thetas = np.empty([numY,numX,2], dtype=np.double)
        thetas[:,:,0], thetas[:,:,1] = np.meshgrid(np.linspace(bottomLeft[0], topRight[0], numX), 
                                                   np.linspace(bottomLeft[1], topRight[1], numY))
        return thetas

    def __init__(self, lensesAndRedshifts, bottomLeft, topRight, numX, numY, renderer = "default", feedbackObject = None, cosmology = "default"):
        """This creates a MultiLensPlane instance that covers the area specified by the `bottomLeft` and
        `topRight` corners. It calculates and stores the deflection angles for the lenses
        in `lensesAndRedshifts`, which should be a list of 
        (:class:`gravitational lens <grale.lenses.GravitationalLens>`, redshift) tuples.
        The deflection angles for the first lens plane are arranges on a grid of `numX` points 
        wide by `numY` points high. The `cosmology` parameter specifies how the redshifts
        should be mapped to angular diameter distances, using a :class:`Cosmology <grale.cosmology.Cosmology>`
        instance.

        If `renderer` is ``None``, the mapping is calculated single threaded, within this
        Python process. Other renderers can be specified as well, for example to calculate the
        mapping faster using multiple cores with the MPI renderer. See the :mod:`grale.renderers`
        module for more information. 
        
        Feedback while rendering can be provided by specifying a `feedbackObject` parameter. See
        the :mod:`grale.feedback` module for more information about allowed values.
        """
        if not lensesAndRedshifts:
            raise MultiLensPlaneException("No lenses and redshifts were specified")

        cosmology = privutil.initCosmology(cosmology)
        if not cosmology:
            raise MultiLensPlaneException("Cosmological model is not set")

        # check Dd vs z/cosmology. This also checks that the lens is probably a lens
        # (i.e. that it has the method getLensDistance())
        count = 0
        for lens, z in lensesAndRedshifts:
            Dd = cosmology.getAngularDiameterDistance(z)
            if abs((Dd - lens.getLensDistance())/Dd) > 1e-6:
                raise MultiLensPlaneException("The specified redshift {} for lens at index {} is not compatible with the angular diameter distance stored in the lens itself".format(z, count))

            count += 1

        renderer, feedbackObject = privutil.initRendererAndFeedback(renderer, feedbackObject, "LENSPLANE")

        thetas = MultiLensPlane._createThetaGrid(bottomLeft, topRight, numX, numY)
        self._thetas = thetas
        self._bottomLeft = bottomLeft[:2]
        self._topRight = topRight[:2]
        self._numX = numX
        self._numY = numY
        self._cosmology = cosmology

        # Create copy of lensesAndRedshifts, and sort on redshift
        lz = sorted(lensesAndRedshifts[:], key = lambda xz: xz[1])
        N = len(lz)

        # Check that no two redshifts are the same
        for i in range(1, N):
            if lz[i-1][1] == lz[i][1]:
                raise MultiLensPlaneException("At least two lenses have the same redshift. Combine them in a CompositeLens instance first.")

        self._caches = _MultiPlaneCache(lz, thetas, renderer, feedbackObject, cosmology, True)
        self._lz = lz

    def getLensesAndRedshifts(self):
        """Returns a copy of the lens and redshift tuples that was specified
        during initialization."""
        return self._lz[:]

    def getRenderInfo(self):
        """Returns a dictionary with the following entries:

         - ``bottomleft``: the bottom-left corner that was specified in the constructor
           of this instance
         - ``topright``: the top-right corner that was specified in the constructor of
           this instance
         - ``xpoints``: the number of points in the x-direction, between the left and right
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled
         - ``ypoints``: the number of points in the y-direction, between the bottom and top
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled
        """
        ri = {
            "bottomleft": self._bottomLeft[:],
            "topright": self._topRight[:],
            "xpoints": self._numX,
            "ypoints": self._numY
        }
        return ri

class MultiImagePlane(object):

    def __init__(self, multiLensPlane, sourceRedshift):
        """Based on the :class:`MultiLensPlane` instance in `multiLensPlane`, which contains the deflection
        fields for a number of lenses at different redshifts, a MultiImagePlane instance
        is created for a certain source plane. This source plane has redshift `sourceRedshift`,
        which will be converted to angular diameter distances using the cosmological
        model stored in `multiLensPlane`.
        """

        # Keep only the part that's relevant to us
        self._zs = sourceRedshift
        self._multiLensPlane = multiLensPlane

        Ti, Tixx, Tixy, Tiyx, Tiyy = multiLensPlane._caches.addSourcePlane(sourceRedshift)

        # Keep only the part that's relevant to us
        lz = multiLensPlane._lz[:]
        while lz and sourceRedshift < lz[-1][1]:
            del lz[-1]
        self._lz = lz

        self._zs = sourceRedshift
        self._betas = Ti
        self._betaderivs = [ Tixx, Tixy, Tiyx, Tiyy ]
        self._invmag = None
        self._invmagGrid = None
        self._critLines = None
        self._betaGrid = None
        self._caustics = None

    def getLensPlane(self):
        return self._multiLensPlane

    def getSourceRedshift(self):
        """Returns the source redshift that was specified during initialization."""
        return self._zs

    def _checkInvMag(self):
        if self._invmag is None:
            betaderivs = self._betaderivs
            self._invmag = betaderivs[0] * betaderivs[3] - betaderivs[1]*betaderivs[2]

    def _calcCrit(self):
        if self._critLines is None:
            self._checkInvMag()
            mlp = self._multiLensPlane
            ctrFinder = contourfinder.ContourFinder(self._invmag, mlp._bottomLeft, mlp._topRight)
            self._critLines = ctrFinder.findContour(0)

    def getCriticalLines(self):
        """This returns a list describing the critical lines associated with this
        image plane. Each entry in the list is itself a list of 2D points, describing
        a connected part of a critical line.
        """
        self._calcCrit()
        return copy.deepcopy(self._critLines)

    def _getBetaGrid(self):
        if not self._betaGrid:
            betas = self._betas
            mlp = self._multiLensPlane
            gx = gridfunction.GridFunction(betas[:,:,0], mlp._bottomLeft, mlp._topRight)
            gy = gridfunction.GridFunction(betas[:,:,1], mlp._bottomLeft, mlp._topRight)
            self._betaGrid = [ gx, gy ]
        else:
            gx, gy = self._betaGrid

        return gx, gy

    def getCaustics(self, approx = False):
        """This returns a list describing the caustics associated with this
        image plane. Each entry in the list is itself a list of 2D points, describing
        a connected part of a caustic.
        """
        self._calcCrit()

        traceFunction = self.traceThetaApproximately if approx else self.traceTheta

        if self._caustics is None:
            gx, gy = self._getBetaGrid()
            contourParts = self._critLines

            self._caustics = [ ]
            for part in contourParts:
                caustPart = traceFunction(part)
                self._caustics.append(np.array(caustPart))

        return copy.deepcopy(self._caustics)

    def traceBetaApproximately(self, beta):
        """Estimates the image plane positions to which the source plane position `beta`
        corresponds. Returns a list of 2D points.
        """
        mlp = self._multiLensPlane
        theta = images.ImagePlane.static_traceBetaApproximately(beta, self._betas, mlp._bottomLeft, mlp._topRight)
        return theta

    def getRenderInfo(self):
        """Returns a dictionary with the following entries:

         - ``bottomleft``: the bottom-left corner that is relevant for this instance. This
           is taken from the LensPlane instance specified in the constructor.
         - ``topright``: the top-right corner that is relevant for this instance. This
           is taken from the LensPlane instance specified in the constructor.
         - ``xpoints``: the number of points in the x-direction, between the left and right
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled. This is taken from the LensPlane
           instance specified in the constructor.
         - ``ypoints``: the number of points in the y-direction, between the bottom and top
           coordinates specified by ``bottomleft`` and ``topright``, at which
           the deflection field was sampled. This is taken from the LensPlane
           instance specified in the constructor.
         - ``xpixels``: based on the image plane to source plane mappings that are known at
           ``xpoints`` * ``ypoints`` grid points, a number of pixels can be defined
           that can contain light from the source plane, and these pixels will
           be used when rendering the image plane or the source plane with
           :func:`renderImages` or :func:`renderSources`. This ``xpixels`` value
           specifies the number of pixels in the x-direction, and is one less
           than ``xpoints``.
         - ``ypixels``: similar to ``xpixels``, but for the y-direction.
        """
        ri = copy.deepcopy(self._multiLensPlane.getRenderInfo())
        ri["xpixels"] = ri["xpoints"]-1
        ri["ypixels"] = ri["ypoints"]-1
        return ri

    def _renderSourcesOrImages(self, sourceList, plane, subSamples, thetaMapping):
        mlp = self._multiLensPlane
        bl, tr = mlp._bottomLeft, mlp._topRight
        numX, numY = mlp._numX, mlp._numY
        
        dstPlane = np.zeros((numY-1, numX-1), dtype=np.double)

        numSub = round(subSamples**0.5)
        if numSub < 1:
            numSub = 1

        dx = np.zeros((numY-1, numX-1, 2), dtype=np.double)
        dy = np.zeros((numY-1, numX-1, 2), dtype=np.double)
        dx[:,:,0] = ((tr[0]-bl[0])/(numX-1))/numSub
        dy[:,:,1] = ((tr[1]-bl[1])/(numY-1))/numSub

        thetas = mlp._thetas[:-1,:-1,:]

        for src in sourceList:
            srcPlane = np.zeros((numY-1, numX-1), dtype=np.double)

            for iy in range(numSub):
                yThetas = thetas + dy*(0.5+iy)
                for ix in range(numSub):
                    xyThetas = yThetas + dx*(0.5+ix) 
                    srcPlane += src.getIntensity(thetaMapping(xyThetas))

            srcPlane /= numSub*numSub
            dstPlane += srcPlane

        if plane is None:
            plane = dstPlane
        else:
            np.copyto(plane, dstPlane)
        return plane

    def _thetaMappingIdentity(self, thetas):
        return thetas

    def traceThetaApproximately(self, thetas):
        """Use the already calculated theta/beta mapping (image plane position to
        source plane positions), to estimate the mapping for theta vectors that
        have not been calculated exactly."""

        betaApprox = np.empty(thetas.shape)
        gx, gy = self._getBetaGrid()

        thetasReshaped = thetas.reshape([-1,2])
        betaApproxReshaped = betaApprox.reshape([-1,2])
        betaApproxReshaped[:,0] = gx.evaluate(thetasReshaped)
        betaApproxReshaped[:,1] = gy.evaluate(thetasReshaped)

        #with open("/tmp/thetalog.txt", "at") as f:
        #    for y in range(thetas.shape[0]):
        #        for x in range(thetas.shape[1]):
        #            f.write("{:g} {:g}\n".format(thetas[y,x,0], thetas[y,x,1]))
        return betaApprox

    def getBetaAndDerivatives(self, thetas):
        """For each theta position in `thetas`, calculate the corresponding position
        in the source plane (beta vector) as well as the derivatives. Note that this is 
        all calculated using a single processor core, no speedup using e.g. 
        OpenMP will be performed."""

        origShape = thetas.shape
        thetas = thetas.reshape((-1,1,2)) # Number of rows differs from one (I think parallelization was done over rows)

        renderer, feedbackObject = privutil.initRendererAndFeedback(None, "none", "LENSPLANE")
        cache = _MultiPlaneCache(self._lz, thetas, renderer, feedbackObject, 
                                 self._multiLensPlane._cosmology, True)

        Ti, Tixx, Tixy, Tiyx, Tiyy = cache.addSourcePlane(self._zs)

        betas = Ti.reshape(origShape)
        Tiall = np.empty((thetas.shape[0], 1, 2, 2))
        Tiall[:,0,0,0] = Tixx[:,0]
        Tiall[:,0,0,1] = Tixy[:,0]
        Tiall[:,0,1,0] = Tiyx[:,0]
        Tiall[:,0,1,1] = Tiyy[:,0]

        if len(origShape) == 1: # Just a single point
            Tiall = Tiall.reshape(2,2)
        else:
            Tiall = Tiall.reshape(origShape[:-1] + (2,2))

        return Ti.reshape(origShape), Tiall
    
    def traceTheta(self, thetas):
        """For each theta position in `thetas`, calculate the corresponding position
        in the source plane. Note that this is all calculated using a single processor
        core, no speedup using e.g. OpenMP will be performed."""

        origShape = thetas.shape
        thetas = thetas.reshape((-1,1,2)) # Number of rows differs from one (I think parallelization was done over rows)

        renderer, feedbackObject = privutil.initRendererAndFeedback(None, "none", "LENSPLANE")
        cache = _MultiPlaneCache(self._lz, thetas, renderer, feedbackObject, 
                                 self._multiLensPlane._cosmology, False)

        Ti = cache.addSourcePlane(self._zs)
        return Ti.reshape(origShape)

    def renderSources(self, sourceList, plane = None, subSamples = 9):
        """For the list of :class:`SourceImage` derived classes in `sourceList`, this function
        calculates what the sources look like based on the dimensions and number of pixels for
        this ImagePlane instance. This is what the image plane would look like if the
        gravitational lens effect could be turned off.

        The function returns a 2D NumPy array containing ``ypixels`` rows, each of ``xpixels``
        pixels wide (see also :func:`getRenderInfo`). If `plane` is specified, the results are
        stored in that 2D NumPy instance, which must have the same dimensions.

        Each pixel is sub-sampled ``sqrt(subSamples)`` times in x- and y- direction, to be able to
        roughly approximate the integration that's needed over the surface area of a pixel.

        Note that a call to only
        `imshow <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.imshow.html>`_.
        will plot the 0,0 value in `plane` (the bottom-left value) as the top-left corner
        causing the result to appear mirrored in the y-direction (the y-axis will point down).
        A subsequent call to `invert_yaxis <https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.invert_yaxis.html>`_
        might be useful.
        """
        return self._renderSourcesOrImages(sourceList, plane, subSamples, self._thetaMappingIdentity)

    def renderImages(self, sourceList, plane = None, subSamples = 9):
        """For the list of :class:`SourceImage` derived classes in `sourceList`, this function
        calculates what the images look like based on the dimensions and number of pixels for
        this ImagePlane instance. 

        The function returns a 2D NumPy array containing ``ypixels`` rows, each of ``xpixels``
        pixels wide (see also :func:`getRenderInfo`). If `plane` is specified, the results are
        stored in that 2D NumPy instance, which must have the same dimensions.

        Each pixel is sub-sampled ``sqrt(subSamples)`` times in x- and y- direction, to be able to
        roughly approximate the integration that's needed over the surface area of a pixel.

        Note that a call to only
        `imshow <https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.imshow.html>`_.
        will plot the 0,0 value in `plane` (the bottom-left value) as the top-left corner
        causing the result to appear mirrored in the y-direction (the y-axis will point down).
        A subsequent call to `invert_yaxis <https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.invert_yaxis.html>`_
        might be useful.
        """
        return self._renderSourcesOrImages(sourceList, plane, subSamples, self.traceThetaApproximately)

    def segment(self, plane, threshold = 0.0):
        """For the image plane `plane` that was rendered using :func:`renderImages`,
        this function looks at all the pixels that have a value larger than `threshold`.
        These pixels are divided into regions that are coherent, and a list of these
        regions is returned. Each region is itself a list of 2D coordinates describing
        the centers of the pixels.
        """
        mlp = self._multiLensPlane
        bl, tr = mlp._bottomLeft, mlp._topRight
        return images.ImagePlane.static_segment(plane, bl, tr, threshold)

    def getInverseMagnificationApproximately(self, thetas):
        """Based on the exact inverse magnifications calculated for theta vectors
        on a grid, calculate the inverse magnifications approximately for arbitrary
        theta-vectors `thetas`."""
        self._checkInvMag()
        if self._invmagGrid is None:
            mlp = self._multiLensPlane
            self._invmagGrid = gridfunction.GridFunction(self._invmag, mlp._bottomLeft, mlp._topRight)

        return self._invmagGrid.evaluate(thetas)

