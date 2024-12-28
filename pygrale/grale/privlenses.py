from .privimages import _getLinesFromInputData
import numpy as np

def _getLenstoolPotentialInfoFromLines(lines):
    potentialInfo = [ ]
    curInfo = { }
    for l in lines:
        if l.startswith("potentiel") or l.startswith("potential"):
            if curInfo:
                potentialInfo.append(curInfo)
            curInfo = { }

        elif l.startswith(" ") or l.startswith("\t"):
            if curInfo is not None:
                parts = l.strip().split()
                if parts[0] == "end":
                    potentialInfo.append(curInfo)
                    curInfo = None
                else:
                    key, value = parts
                    curInfo[key] = float(value)
        else:
            if curInfo:
                potentialInfo.append(curInfo)
            curInfo = None

    if curInfo:
        potentialInfo.append(curInfo)
    return potentialInfo

def _getLenstoolCosmologyFromLines(lines):
    from .lenses import LensException

    idx = None
    for i in range(len(lines)):
        if lines[i].startswith("cosmologie") or lines[i].startswith("cosmology"):
            idx = i
            break
    else:
        raise LensException("No cosmology settings found in file")

    settings = { "omegaK": 0 } # Doesn't always seem to be present
    while True:
        idx += 1
        l = lines[idx]
        if not (l.startswith(" ") or l.startswith("\t")):
            break
        parts = l.strip().split()
        if parts[0] == "end":
            break
        settings[parts[0]] = float(parts[1])

    from . import cosmology
    return cosmology.Cosmology(settings["H0"]/100.0,
                               settings["omegaM"],
                               settings["omegaK"],
                               settings["omegaX"],
                               settings["wX"])

def createLensFromLenstoolFile(inputData, mirrorX = False, cosmology = None):
    """Based on a `LensTool <https://projets.lam.fr/projects/lenstool/wiki>`_ model,
    a corresponding :class:`lens model<grale.lenses.GravitationalLens>` is constructed.
    The function returns a tuple consisting of the lens model, the lens redshift and
    the cosmological model. An illustration can be found in the notebook
    `lenstooltest.ipynb <_static/lenstooltest.ipynb>`_

    Note that this is preliminary code, and currently only the PIEMD model (LensTool
    model type 81) is handled.

    Arguments:

     - `inputData`: the data from a LensTool file (typically with '.par' extension), either
       as a file name, a file object or a string.
     - `mirrorX`: if ``True``, the model will be mirrored along the X-axis.
     - `cosmology`: if specified, the cosmological model from the input file will be ignored
       and this one will be used.
    """
    lines = _getLinesFromInputData(inputData)
    potentialInfo = _getLenstoolPotentialInfoFromLines(lines)
    cosm = _getLenstoolCosmologyFromLines(lines) if not cosmology else cosmology
    if not potentialInfo:
        return (None, None, cosm)

    from .constants import ANGLE_ARCSEC, CONST_G, DIST_KPC
    from . import lenses

    LensException = lenses.LensException

    def _ltPIEMDHandler(info, Dd):
        epsHat = info["ellipticite"] if "ellipticite" in info else info["ellipticity"]
        a = info["core_radius"]*ANGLE_ARCSEC if "core_radius" in info else info["core_radius_kpc"]/(Dd/DIST_KPC)
        s = info["cut_radius"]*ANGLE_ARCSEC if "cut_radius" in info else info["cut_radius_kpc"]/(Dd/DIST_KPC)
        sigma = info["v_disp"]*1000

        centralDensity = lenses.PIEMDLens.getCentralDensityFromVelocityDispersion(sigma, a, s, Dd)
        eps = lenses.PIEMDLens.getEpsilonFromEllipticity(epsHat)

        if eps == 0: # Can't use PIEMDLens for this
            return lenses.PIMDLens(Dd, { "centraldensity": centralDensity, "coreradius": a, "scaleradius": s })

        return lenses.PIEMDLens(Dd, { "centraldensity": centralDensity, "coreradius": a, "scaleradius": s,
                                      "epsilon": eps })

    _lenstoolPotentialHandlers = { 81: _ltPIEMDHandler }

    subLenses = [ ]
    for p in potentialInfo:
        profileId = int(p["profil"]) if "profil" in p else p["profile"]
        if not profileId in _lenstoolPotentialHandlers:
            raise LensException("No handler for profile {}".format(profileId))

        handler = _lenstoolPotentialHandlers[profileId]
        zd = p["z_lens"]
        Dd = cosm.getAngularDiameterDistance(zd)
        lens = handler(p, Dd)

        x = p["x_centre"]*ANGLE_ARCSEC
        y = p["y_centre"]*ANGLE_ARCSEC
        angle = p["angle_pos"] # TODO: check!
        
        if mirrorX:
            x = -x
            angle = 180-angle

        subLenses.append({ "lens": lens, "x": x, "y": y, "factor": 1, "angle": angle})

    Dd_first = subLenses[0]["lens"].getLensDistance()
    for d in subLenses:
        if d["lens"].getLensDistance() != Dd_first:
            raise LensException("Not all lens components have same lens distance")

    return (lenses.CompositeLens(Dd_first, subLenses), zd, cosm)

def _getFactorListFromKernel(K, scale, diOffset, djOffset):
    res = []
    for y in range(K.shape[0]):
        i = -K.shape[0]//2 + 1 + y
        for x in range(K.shape[1]):
            j = -K.shape[1]//2 + 1 + x

            val = -K[y,x]*scale # Changed a sign convention somewhere, hence the '-'; TODO: fix this
            if val != 0:
                res.append({ "factor": val, "di": i+diOffset, "dj": j+djOffset })

    return res

def _simplifyFactorList(l):
    factors = { }
    for d in l:
        didj = d["di"],d["dj"]
        if didj in factors:
            factors[didj] += d["factor"]
        else:
            factors[didj] = d["factor"]

    r = [ { "factor": factors[didj], "di": didj[0], "dj": didj[1] } for didj in factors if factors[didj] != 0 ]
    #print("Simplified:")
    #import pprint
    #pprint.pprint(r)
    return r

def _getDensityScaleFactor(thetas, Dd, phiScale, laplacianKernel):
    from . import lenses
    # Calculate the factor needed to convert density according to the laplacian
    # TODO: no need to use all thetas
    phiSheet = lenses.MassSheetLens(Dd, { "density": 1.0 }).getProjectedPotential(1,1,thetas)
    phiSheet -= np.min(phiSheet)
    phiSheet /= phiScale # Use same scale factor
    # Convolve just a small part
    return np.sum(phiSheet[:laplacianKernel.shape[0],:laplacianKernel.shape[1]] * laplacianKernel)  

def _getDeflectionScaleFactor(thetas, Dd, phiScale):
    from . import lenses
    from .constants import SPEED_C
    # A SIS lens has a contant deflection
    velDisp = 400000
    sisLens = lenses.SISLens(Dd, { "velocityDispersion": velDisp})
    sisPot = sisLens.getProjectedPotential(1,1,thetas)
    sisPot /= phiScale

    expectedDeflection = 4*np.pi*velDisp**2/SPEED_C**2 # For Dds/Ds = 1
    
    potGradKernelX = np.array([[1.0, -1.0]])
    potGradKernelY = np.array([[1.0], [-1.0]])
    
    import scipy.ndimage
    sisAx = scipy.ndimage.convolve(sisPot, potGradKernelX)[1:-1,1:-1]
    sisAy = scipy.ndimage.convolve(sisPot, potGradKernelY)[1:-1,1:-1]
    # TODO: this is not entirely correct I think, since the Ax and Ay components don't
    #       really refer to the same point. Seems good enough though
    sisAsize = (sisAx**2 + sisAy**2)**0.5
    goodPos = ~np.isnan(sisAsize)
    sisAstd = np.std(sisAsize[goodPos])
    sisAsize = np.mean(sisAsize[goodPos])
    
    return sisAsize/expectedDeflection # need to multiply real alphas with this to get convolution results

def createEquivalentPotentialGridLens(lens, bottomLeft, topRight, NX, NY, maskRegions,
                                      potentialGradientWeight, densityGradientWeight,
                                      densityWeight, pixelEnlargements=2,
                                      enlargeDiagonally=False, circleToPolygonPoints=10000,
                                      feedbackObject="default", qpsolver="scs",
                                      laplacianKernel = np.array([[  0, 0,  1, 0, 0 ],
                                                                  [  0, 1,  2, 1, 0 ],
                                                                  [  1, 2,-16, 2, 1 ],
                                                                  [  0, 1,  2, 1, 0 ],
                                                                  [  0, 0,  1, 0, 0 ]],dtype=np.double),
                                      # No pixel enlargements are used for these masks!
                                      #  [ { "maskRegions": ..., "density": ..., ("upperlimit": ...) }, ... ]
                                      maxDensityConstraints = [],
                                      #  [ { "maskRegions": ..., "density": ...}, ... ]
                                      exactDensityConstraints = [],
                                      # [ { "maskRegions": ..., "ax": ..., "ay": ... }] # Here you need to take into account that gradients are actually between pixels
                                      # or
                                      # [ { "maskRegions": ..., "lens": ... }]
                                      exactDeflectionConstraints = [],
                                      exactDeflectionTolerance = 0,
                                      ignorePixelMismatch = False,
                                      ignorePositiveDensityConstraint = False,
                                      exceptionOnFail = True,
                                      dryrun = False
                                      ):
    r"""This uses a quadratic programming approach to extrapolate the lens potential values
    in certain regions (typically covering the images in a lensing system), thereby creating
    a lens that has the same effect (because the lens potential is the same in the image
    regions).

    Examples can be found in the `msdexample-equivlenstests.ipynb <_static/msdexample-equivlenstests.ipynb>`_
    and `potentialextrap_multisheet.ipynb <_static/potentialextrap_multisheet.ipynb>`_ notebooks.

    Arguments:
     - `lens`: the procedure will start from the lens potential values of this
       :class:`GravitationalLens <grale.lenses.GravitationalLens>` model. It will
       sample the lens potential on a grid, and keep some values fixed, determined
       by a mask. The other lens potential values of the new lens model will be
       extrapolated.
     - `bottomLeft`, `topRight`, `NX`, `NY`: these values determine the grid
       on which the lens potential will be sampled from the original lens, and
       which will be used to define the new lens (a :class:`PotentialGridLens <grale.lenses.PotentialGridLens>`).
       See :func:`createThetaGrid <grale.util.createThetaGrid>`.
     - `maskRegions`: this is a list of region descriptions that will be combined
       to create a binary mask of `NX` by `NY` that indicates which lens potential
       values should be kept. This is passed to the :func:`createThetaGridAndImagesMask <grale.util.createThetaGridAndImagesMask>`
       function.
     - `potentialGradientWeight`, `densityGradientWeight`, `densityWeight`: the
       quadratic programming problem tries to optimize a combination of three parts:
       one for the gradient of the lens potential, one for the gradient of the resulting
       mass density, and one for the mass density itself. These weights are used to specify
       their respective contributions.
     - `pixelEnlargements`, `enlargeDiagonally`, `circleToPolygonPoints`: these are
       passed on to the similarly named arguments of :func:`createThetaGridAndImagesMask <grale.util.createThetaGridAndImagesMask>`.
     - `feedbackObject`: can be used to specify a particular :ref:`feedback mechanism <feedback>`.
     - `qpsolver`: for the quadratic programming optimization, the `qpsolvers <https://pypi.org/project/qpsolvers/>`_
       module is used, which itself allows for different solver implementations to be used.
       This argument specifies the name of the solver, and is passed as the `solver`
       argument to `qpsolvers.solve_qp <https://qpsolvers.github.io/qpsolvers/quadratic-programming.html#qpsolvers.solve_qp>`_.
     - `laplacianKernel`: to go from the lens potential values to the density, a convolution
       with this 2D kernel is performed. The default is

        .. math::

            \left[
                \begin{array}{cccc}
                    0 & 0 & 1 & 0 & 0 \\
                    0 & 1 & 2 & 1 & 0 \\
                    1 & 2 & -16 & 2 & 1 \\
                    0 & 1 & 2 & 1 & 0 \\
                    0 & 0 & 1 & 0 & 0 
                \end{array}
            \right]

     - `maxDensityConstraints`: since this kernel allows the density to be calculated at the
       grid points, constraints can be specified for those. This is a list of these constraints
       where each entry is a dictionary with the following keys:

        - `maskRegions`: a list of regions to which the specified density applies. Is passed
          on to the :func:`createThetaGridAndImagesMask <grale.util.createThetaGridAndImagesMask>`
          again, without any pixel enlargements.
        - `density`: can either be a scalar value or a grid with the same dimensions as the
          lens potential grid.
        - `upperlimit` (optional): if omitted, this is interpreted as ``True``, and means that
          the specified density in the specified region should not be exceeded. Set this to
          ``False`` for a lower bound on the density instead.

     - `exactDensityConstraints`: similar to `maxDensityConstraints`, but requests the exact
       mass density values in certain regions. This is again a list of similar dictionaries, but
       in this case only with `maskRegions` and `density` keys.
     - `exactDeflectionConstraints`: instead of providing constraints on the density, constraints
       on the deflection field (the gradient of the lens potential) can be set as well. This is
       a list of dictionaries with entries

         - `maskRegions`: again a list of regions, this time for which the specified deflection
           field applies.
         - `lens`: get the desired deflection field at the grid points from this lens.

        Alternatively, instead of the `lens` key, `ax` and `ay` keys may be present, each describing
        one component of the deflection field. In case these are specified manually, it will be
        necessary to take into account that the kernel used (e.g. [-1, 1]) gives an approximation
        in between grid points.
     - `exactDeflectionTolerance`: using the constraints above as exact constraints rarely works.
       Instead it is possible to request that this deflection field is obtained within some
       tolerance, e.g. 0.1 arcsec, on each side of the exact value.
     - `ignorePixelMismatch`: if the specified grid parameters does not yield the same spacing
       between grid points in x- and y-direction, the routine will abort by default. You can override
       this and continue anyway by setting this flag.
     - `ignorePositiveDensityConstraint`: by default, a constraint is added that requires the
       density resulting from the lens potential to be positive everywhere. To ignore this condition,
       set this flag.
     - `exceptionOnFail`: by default, if the quadratic programming solver is unable to come up with
       a solution, an exception is thrown. Setting this to ``False`` disables this and returns the
       dictionary mentioned below anyway, but with the resulting solution and resulting lens set
       to ``None``.
     - `dryrun`: only build masks and matrices, don't actually try to solve anything.

    The return value is a dictionary with the following keys:

     - `philens_equiv`: this the the main result of the optimization routine, a 
       :class:`PotentialGridLens <grale.lenses.PotentialGridLens>` that's based on lens potential
       values defined on the grid that you specified, where some values were taken from the input
       lens, and others were optimized according to the weights and constraints that were
       specified.
     - `philens_orig`: as an intermediate step, the input lens is approximated by sampling
       the lens potential on the grid points, and this is the approximate lens model. You can
       compare this to the original model to see how much sense the approximation makes.
     - `mask`: the main `maskRegions` argument, together with other parameters (e.g.
       `pixelEnlargements`) leads to this binary mask. Where it's ``True``, the lens potential
       values are fixed, the others will be optimized.
     - `masksMaxDens`: each entry of `maxDensityConstraints` has its own `maskRegions`
       argument, leading to a binary mask again. This entry contains a list of these masks.
     - `masksExactDens`: similar, but for the `exactDensityConstraints`.
     - `masksDeflection`: similar, but for the `exactDeflectionConstraints`.
     - `P`, `q`, `G`, `h`, `A`, `b`: these are the matrices that are created for the
       quadratic programming problem
     - `initvals`: the initial values that are passed to the solver, these are the
       lens potential values of the input lens.
     - `x`: this contains the solution to the quadratic programming problem.

    """
    import time
    from . import feedback
    from . import util
    from . import lenses
    from . import quadprogmatrix
    from . import privutil
    from .constants import ANGLE_ARCSEC
    from .lenses import LensException
    from qpsolvers import solve_qp 
    import scipy.sparse as sparse

    if not ignorePixelMismatch:
        distScale = ((topRight[0]-bottomLeft[0])**2 + (topRight[1]-bottomLeft[1])**2)**0.5
        Nscale = (NX**2 + NY**2)**0.5
        dx = ((topRight[0] - bottomLeft[0])/distScale)*(Nscale/NX)
        dy = ((topRight[1] - bottomLeft[1])/distScale)*(Nscale/NX)
        if abs(dx-dy) > 1e-5:
            raise LensException("Pixel sizes in x- and y-directions don't seem to match well enough, use 'ignorePixelMismatch' to override")

    if NY < laplacianKernel.shape[0] or NX < laplacianKernel.shape[1]:
        raise LensException("Grid is not large enough to use Laplacian kernel")

    if laplacianKernel.shape[0] != laplacianKernel.shape[1] or laplacianKernel.shape[0] == 0 or laplacianKernel.shape[0]%2 != 1:
        raise LensException("Laplacian kernel must have same number of rows and colums, which must be odd")

    t0 = time.time()
    kHw = laplacianKernel.shape[0]//2
    feedbackObject = privutil.processFeedbackObjectArgument(feedbackObject)

    lens2, maskRegions2 = None, None
    if type(lens) is tuple or type(lens) is list:
        if len(lens) != 2:
            raise LensException("For constraints from different lenses, only two are allowed")

        lens, lens2 = lens
        # mask regions should also contain masks for the two lenses
        maskRegions, maskRegions2 = maskRegions

    # The first lens gets us started either way
    thetas, mask = util.createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, maskRegions, pixelEnlargements,
                                                     enlargeDiagonally, circleToPolygonPoints)

    mask = mask.astype(np.intc)
    
    feedbackObject.onStatus("Calculating lens potential values")
    phi = lens.getProjectedPotential(1,1,thetas)
    phi -= np.min(phi)
    phiScale = np.max(phi)

    # If there are two lenses, modify the mask based on the second one. This does not
    # affect the scale factors anymore, we'll just keep those 
    if lens2:
        _, mask2 = util.createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, maskRegions2, pixelEnlargements,
                                                     enlargeDiagonally, circleToPolygonPoints)
        
        # Here, we use 2 for the regions that stay fixed, up to a constant
        # and a gradient
        mask2 = mask2.astype(np.intc) * 2
        mask += mask2

        if np.max(mask) > 2: # Then we've added a 1 to a 2
            raise LensException("Regions for the two lenses seem to overlap, can't construct a total mask")

        # Setup the correct initial potential values for this lens as well
        phi2 = lens2.getProjectedPotential(1,1,thetas)
        phi2 -= np.min(phi)

        phi[mask == 2] = phi2[mask == 2]
        phiScale = max(phiScale, np.max(phi))

    if phiScale == 0:
        phiScale = ANGLE_ARCSEC*ANGLE_ARCSEC

    unitDensityScaleFactor = _getDensityScaleFactor(thetas, lens.getLensDistance(), phiScale, laplacianKernel)
    deflectionAngleScaleFactor = _getDeflectionScaleFactor(thetas, lens.getLensDistance(), phiScale)

    #import pickle
    #pickle.dump(phi, open("debugphi.dat", "wb"))
    prob = quadprogmatrix.MaskedPotentialValues(phi, mask, phiScale)
        
    feedbackObject.onStatus("Calculating linear constraints")
    
    laplacian = _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1.0, 0, 0))
    gradientDi = [ { "factor": 1, "di": 0, "dj": 0}, { "factor": -1, "di": 1, "dj": 0} ]
    gradientDj = [ { "factor": 1, "di": 0, "dj": 0}, { "factor": -1, "di": 0, "dj": 1} ]

    if ignorePositiveDensityConstraint:
        G, h = None, None
    else:
        G, h = prob.getLinearConstraintMatrices(laplacian)
    
    maxMasks, exactMasks, exactGradMasks = [], [], []
    for constr in maxDensityConstraints:
        _, constrMask = util.createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, constr["maskRegions"], 0)
        G2, h2 = prob.getLinearConstraintMatrices(laplacian, relevantGridPositions=constrMask, 
                                                  # The construction below makes sure that there's a 2D array, even if just a single
                                                  # value for the density was specified
                                                  limitingValues=np.ones((NY,NX), dtype=np.double) * constr["density"] * unitDensityScaleFactor,
                                                  isUpperLimit=constr["upperlimit"] if "upperlimit" in constr else True)
        if G is None:
            G, h = G2, h2
        else:
            G = sparse.vstack([G,G2])
            h = np.concatenate([h,h2])
        
        maxMasks.append(constrMask)

    A, b = None, None

    # TODO: for testing, disable gradient by using constraint
    #A = np.zeros((3, prob.getNumberOfVariables()), dtype=np.double)
    #b = np.zeros((3,))
    #A[0,1] = 1
    #A[1,2] = 1
    #A[2,0] = 1
    
    for constr in exactDensityConstraints:
        _, constrMask = util.createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, constr["maskRegions"], 0)
        A2, b2 = prob.getLinearConstraintMatrices(laplacian, relevantGridPositions=constrMask,
                                                  # The construction below makes sure that there's a 2D array, even if just a single
                                                  # value for the density was specified
                                                  limitingValues=np.ones((NY,NX), dtype=np.double) * constr["density"] * unitDensityScaleFactor)
        if A is None:
            A, b = A2, b2
        else:
            A = sparse.vstack([A,A2])
            b = np.concatenate([b,b2])
            
        exactMasks.append(constrMask)
   
    for constr in exactDeflectionConstraints:
        _, constrMask = util.createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, constr["maskRegions"], 0)
        
        if "lens" in constr: # Calculate the deflection angles from this lens
            constrLens = constr["lens"]
            constr = { }
            dxy = (topRight - bottomLeft)/np.array([NX,NY], dtype=np.double)

            # Make sure the measurements are actually _in between_ the pixels
            constr["ax"] = constrLens.getAlphaVector(thetas + np.array([dxy[0], 0],dtype=np.double)/2)[:,:,0]
            constr["ay"] = constrLens.getAlphaVector(thetas + np.array([0, dxy[1]],dtype=np.double)/2)[:,:,1]
            
        for gradientKernel, constrName in [ (gradientDj, "ax"), (gradientDi, "ay")]:
            if constr[constrName] is None: # skip this
                continue

            if not exactDeflectionTolerance: # exact

                A2, b2 = prob.getLinearConstraintMatrices(gradientKernel, relevantGridPositions=constrMask,
                                                          limitingValues=np.ones((NY,NX), dtype=np.double) * constr[constrName] * deflectionAngleScaleFactor)
                if A is None:
                    A, b = A2, b2
                else:
                    A = sparse.vstack([A,A2])
                    b = np.concatenate([b,b2])
        
            else: # allow some tolerance
                
                for upper, tol in [ (True, exactDeflectionTolerance), (False, -exactDeflectionTolerance)]:

                    G2, h2 = prob.getLinearConstraintMatrices(gradientKernel, relevantGridPositions=constrMask, 
                                                      limitingValues=(np.ones((NY,NX), dtype=np.double) * constr[constrName] + tol)* deflectionAngleScaleFactor,
                                                      isUpperLimit=upper)
                    if G is None:
                        G, h = G2, h2
                    else:
                        G = sparse.vstack([G,G2])
                        h = np.concatenate([h,h2])

        exactGradMasks.append(constrMask)
            
    w1 = potentialGradientWeight
    w2 = densityGradientWeight
    w3 = densityWeight

    feedbackObject.onStatus("Calculating quadratic optimization matrices")

    laplacianGradientDi = _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1.0, 0, 0) + _getFactorListFromKernel(laplacianKernel, -1.0, 1, 0))
    laplacianGradientDj = _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1.0, 0, 0) + _getFactorListFromKernel(laplacianKernel, -1.0, 0, 1))
    #laplacianGradientDij = _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1/2**0.5, 0, 0) + _getFactorListFromKernel(laplacianKernel, -1/2**0.5, 1, 1))
    #laplacianGradientDij2 = _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1/2**0.5, 0, 0) + _getFactorListFromKernel(laplacianKernel, -1/2**0.5, 1, -1))

    P, q = None, None
    for weight, kernel in [
            (w1, gradientDi),
            (w1, gradientDj),
            (w2, laplacianGradientDi),
            (w2, laplacianGradientDj),
            #(w2, laplacianGradientDij),
            #(w2, laplacianGradientDij2),
            (w3, laplacian)
        ]:
        if weight == 0:
            continue

        P2, q2 = prob.getQuadraticMinimizationMatrices(kernel)
        P2 *= weight
        q2 *= weight

        if P is None:
            P, q = P2, q2
        else:
            P += P2
            q += q2

    useWarmStart = True
    if qpsolver.lower() == "mosek": # Mosek does not use warm start, this gets rid of warning message
        useWarmStart = False
    initVals = prob.getInitialValues() if useWarmStart else None
    
    feedbackObject.onStatus("Solving quadratic programming problem")

    newLens = None
    sol = None

    if not dryrun:
        if P is None:
            raise LensException("Nothing optimize (all weights zero?)")

        sol = solve_qp(P, q, G, h, A, b, solver=qpsolver, initvals=initVals)

    if sol is None:
        if dryrun:
            pass # Ok, didn't request a solution
        else:
            if exceptionOnFail:
                raise LensException("Unable to solve quadratic programming problem")
            print("WARNING: Unable to solve quadratic programming problem")
    else:
        newPhi = prob.getFullSolution(sol)
        newLens = lenses.PotentialGridLens(lens.getLensDistance(), { "values": newPhi, "bottomleft": bottomLeft, "topright": topRight})
    
    t1 = time.time()
    feedbackObject.onStatus("Done, in {:.3g} seconds".format(t1-t0))

    return {
        "mask": mask,
        "masksMaxDens": maxMasks,
        "masksExactDens": exactMasks,
        "masksDeflection": exactGradMasks,
        "philens_orig": lenses.PotentialGridLens(lens.getLensDistance(), { "values": phi, "bottomleft": bottomLeft, "topright": topRight}),
        "philens_equiv": newLens,
        # QP parameters
        "P": P,
        "q": q,
        "G": G,
        "h": h,
        "A": A,
        "b": b,
        "x": sol,
        "initvals": initVals,
        # scale factors
        "densityscale": unitDensityScaleFactor,
        "deflectionscale": deflectionAngleScaleFactor
    }
