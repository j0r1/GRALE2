from .privimages import _getLinesFromInputData
import numpy as np

def _getLenstoolPotentialInfoFromLines(lines):
    potentialInfo = [ ]
    curInfo = { }
    for l in lines:
        if l.startswith("potentiel"):
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
        if lines[i].startswith("cosmologie"):
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

def createLensFromLenstoolFile(inputData, mirrorX = False):
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
    """
    lines = _getLinesFromInputData(inputData)
    potentialInfo = _getLenstoolPotentialInfoFromLines(lines)
    cosm = _getLenstoolCosmologyFromLines(lines)
    if not potentialInfo:
        return (None, None, cosm)

    from .constants import ANGLE_ARCSEC, CONST_G, DIST_KPC
    from . import lenses

    LensException = lenses.LensException

    def _ltPIEMDHandler(info, Dd):
        epsHat = info["ellipticite"]
        a = info["core_radius"]*ANGLE_ARCSEC if "core_radius" in info else info["core_radius_kpc"]/(Dd/DIST_KPC)
        s = info["cut_radius"]*ANGLE_ARCSEC if "cut_radius" in info else info["cut_radius_kpc"]/(Dd/DIST_KPC)
        sigma = info["v_disp"]*1000

        centralDensity = (3*sigma**2)/(4*CONST_G*Dd) * (s**2-a**2)/(a*s**2)
        e = 1.0-((1.0-epsHat)/(1.0+epsHat))**0.5
        eps = e/(2-e)
        #print("centralDensity", centralDensity)
        #print("e", e)
        #print()

        return lenses.PIEMDLens(Dd, { "centraldensity": centralDensity, "coreradius": a, "scaleradius": s,
                                      "epsilon": eps })

    _lenstoolPotentialHandlers = { 81: _ltPIEMDHandler }

    subLenses = [ ]
    for p in potentialInfo:
        profileId = int(p["profil"])
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
                                      ignorePixelMismatch = False,
                                      ignorePositiveDensityConstraint = False,
                                      exactDeflectionTolerance = 0,
                                      exceptionOnFail = True
                                      ):

    """This uses a quadratic programming approach to extrapolate the lens potential values
    in certain regions (typically covering the images in a lensing system), thereby creating
    a lens that has the same effect (because the lens potential is the same in the image
    regions).

    Arguments:
     - `lens`
     - `bottomLeft`, `topRight`, `NX`, `NY`
     - `maskRegions`
     - `potentialGradientWeight`, `densityGradientWeight`, `densityWeight`
     - `pixelEnlargements`, `enlargeDiagonally`, `circleToPolygonPoints`
     - `feedbackObject`
     - `qpsolver`
     - `laplacianKernel`
     - `maxDensityConstraints`
     - `exactDensityConstraints`
     - `exactDeflectionConstraints`
     - `ignorePixelMismatch`
     - `ignorePositiveDensityConstraint`
     - `exactDeflectionTolerance`
     - `exceptionOnFail`
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

    thetas, mask = util.createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, maskRegions, pixelEnlargements,
                                                     enlargeDiagonally, circleToPolygonPoints)
    
    feedbackObject.onStatus("Calculating lens potential values")
    phi = lens.getProjectedPotential(1,1,thetas)
    phi -= np.min(phi)
    phiScale = np.max(phi)
    
    unitDensityScaleFactor = _getDensityScaleFactor(thetas, lens.getLensDistance(), phiScale, laplacianKernel)
    deflectionAngleScaleFactor = _getDeflectionScaleFactor(thetas, lens.getLensDistance(), phiScale)
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

    P, q = None, None
    for weight, kernel in [
            (w1, gradientDi),
            (w1, gradientDj),
            (w2, laplacianGradientDi),
            (w2, laplacianGradientDj),
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

    if P is None:
        raise LensException("Nothing optimize (all weights zero?)")

    useWarmStart = True
    if qpsolver.lower() == "mosek": # Mosek does not use warm start, this gets rid of warning message
        useWarmStart = False
    initVals = prob.getInitialValues() if useWarmStart else None
    
    feedbackObject.onStatus("Solving quadratic programming problem")

    newLens = None
    sol = solve_qp(P, q, G, h, A, b, solver=qpsolver, initvals=initVals) # Hard constraints
    if sol is None:
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
    }
