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

            val = -K[y,x]*scale # Changed a sign convention somewhere, hence the '-'
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
                                      ignorePixelMismatch = False
                                      ):

    """TODO"""
    import time
    from . import feedback
    from . import util
    from . import lenses
    from . import quadprogmatrix
    from .constants import ANGLE_ARCSEC
    from .lenses import LensException
    from qpsolvers import solve_qp 

    if not ignorePixelMismatch:
        distScale = ((topRight[0]-bottomLeft[0])**2 + (topRight[1]-bottomLeft[1])**2)**0.5
        Nscale = (NX**2 + NY**2)**0.5
        dx = ((topRight[0] - bottomLeft[0])/distScale)*(Nscale/NX)
        dy = ((topRight[1] - bottomLeft[1])/distScale)*(Nscale/NX)
        if abs(dx-dy) > 1e-5:
            raise LensException("Pixel sizes in x- and y-directions don't seem to match well enough, use 'ignorePixelMismatch' to override")

    t0 = time.time()

    if feedbackObject is None:
        feedbackObject = feedback.Feedback()
    else:
        if type(feedbackObject) == str:
            feedbackObject = feedback.getFeedbackClass(feedbackObject)
            feedbackObject = feedbackObject()
        else:
            # Assume this is a created instance
            pass

    thetas, mask = util.createThetaGridAndImagesMask(bottomLeft, topRight, NX, NY, maskRegions, pixelEnlargements,
                                                     enlargeDiagonally, circleToPolygonPoints)

    feedbackObject.onStatus("Calculating lens potential values")
    phi = lens.getProjectedPotential(1,1,thetas)
    phiLens = lenses.PotentialGridLens(lens.getLensDistance(), { "values": phi, "bottomleft": bottomLeft, "topright": topRight})

    phiMin = np.min(phi)
    phi -= phiMin
    phiScale = np.max(phi)
    
    prob = quadprogmatrix.MaskedPotentialValues(phi, mask, phiScale)
    
    feedbackObject.onStatus("Calculating linear constraints")
    G,h = prob.getLinearConstraintMatrices(_simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1.0, 0, 0)))
    
    w1 = potentialGradientWeight
    w2 = densityGradientWeight
    w3 = densityWeight

    feedbackObject.onStatus("Calculating quadratic optimization matrices")
    P, q = prob.getQuadraticMinimizationMatrices([
        { "weight": w1, "kernel": [
            { "factor": 1, "di": 0, "dj": 0}, { "factor": -1, "di": 0, "dj": 1}
        ] },
        { "weight": w1, "kernel": [
            { "factor": 1, "di": 0, "dj": 0}, { "factor": -1, "di": 1, "dj": 0}
        ] },
        { "weight": w2, "kernel": _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1.0, 0, 0) + _getFactorListFromKernel(laplacianKernel, -1.0, 0, 1)) },
        { "weight": w2, "kernel": _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1.0, 0, 0) + _getFactorListFromKernel(laplacianKernel, -1.0, 1, 0)) },
        { "weight": w3, "kernel": _simplifyFactorList(_getFactorListFromKernel(laplacianKernel, 1.0, 0, 0)) },
    ])

    initVals = prob.getInitialValues()
    
    feedbackObject.onStatus("Solving quadratic programming problem")
    sol = solve_qp(P, q, G, h, solver=qpsolver, initvals=initVals)
    if sol is None:
        raise LensException("Unable to solve quadratic programming problem")
    
    newPhi = prob.getFullSolution(sol)
    newPhiLens = lenses.PotentialGridLens(lens.getLensDistance(), { "values": newPhi, "bottomleft": bottomLeft, "topright": topRight})
    
    t1 = time.time()
    feedbackObject.onStatus("Done, in {:.3g} seconds".format(t1-t0))

    return {
        "mask": mask,
        "philens_orig": phiLens,
        "philens_equiv": newPhiLens,
        # QP parameters
        "P": P,
        "q": q,
        "G": G,
        "h": h,
        "x": sol,
        "initvals": initVals,
    }

