from .privimages import _getLinesFromInputData

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
    idx = None
    for i in range(len(lines)):
        if lines[i].startswith("cosmologie"):
            idx = i
            break
    else:
        raise Exception("No cosmology settings found in file")

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
    lines = _getLinesFromInputData(inputData)
    potentialInfo = _getLenstoolPotentialInfoFromLines(lines)
    cosm = _getLenstoolCosmologyFromLines(lines)
    if not potentialInfo:
        return (None, None, cosm)

    from .constants import ANGLE_ARCSEC, CONST_G, DIST_KPC
    from . import lenses

    def _ltPIEMDHandler(info, Dd):
        eps = info["ellipticite"]
        a = info["core_radius"]*ANGLE_ARCSEC if "core_radius" in info else info["core_radius_kpc"]/(Dd/DIST_KPC)
        s = info["cut_radius"]*ANGLE_ARCSEC if "cut_radius" in info else info["cut_radius_kpc"]/(Dd/DIST_KPC)
        sigma = info["v_disp"]*1000

        centralDensity = (3*sigma**2)/(4*CONST_G*Dd) * (s**2-a**2)/(a*s**2)
        #eps = 1.0-((1.0-eps)/(1.0+eps))**0.5 # TODO: is this necessary?
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
            raise Exception("No handler for profile {}".format(profileId))

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
            raise Exception("Not all lens components have same lens distance")

    return (lenses.CompositeLens(Dd_first, subLenses), zd, cosm)

