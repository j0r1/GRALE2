"""This is a module with a helper function to import parametric inversion
settings from a Lenstool configuration file. The module itself can also
be executed on the command line to generate template inversion scripts
based on a Lenstool configuration.

Usage::

    python -m grale.ltparamdesc ...

"""
from .constants import *
import sys
import os

# "x", "y", "angle", "velocitydispersion", "coreradius", "scaleradius", "ellipticity" entries
def _addBlock(outputLines, block, factor=1):
    try:
        ellStr = block["ellipticity"].strip()
        if ellStr[-1] == ",":
            ellStr = ellStr[:-1]
        ellipt = float(ellStr)
        if ellipt == 0:
            # TODO: create PIMDLens/LTPIMDLens OpenCL version and use this
            block["ellipticity"] = f"0.0001, # Was really '{block['ellipticity']}' but can't use this in PIEMD lens"
    except Exception as e:
        #print(e, file=sys.stderr)
        pass

    comment = ""
    if block['comment']:
        comment = "# " + block['comment'] + "\n        "
    o = f"""    {{
        {comment}"x": {block['x']}
        "y": {block['y']}
        "angle": {block['angle']}
        "factor": {factor},
        "lens": {{
            "type": "LTPIEMDLens",
            "params": {{
                "velocitydispersion": {block['velocitydispersion']}
                "coreradius": {block['coreradius']}
                "scaleradius": {block['scaleradius']}
                "ellipticity": {block['ellipticity']}
            }}
        }}
    }},"""
    for l in o.splitlines():
        outputLines.append(l)

def _processPotfileFile(fileName, mag0, slope, vdSlope, coreStr, bsSigma, bsCut, vdVarName, cutVarName,
                        useRADirection, lensDistanceString, arcsecString, degreeString, kpcString):

    from . import images

    outputLines = []

    lines = open(fileName, "rt").readlines()
    firstLine = lines[0]
    del lines[0]

    parts = firstLine.strip().split()
    if parts[0] != "#REFERENCE":
        raise Exception(f"Expecting galaxy file '{fileName}' to start with '#REFERENCE'")

    if int(parts[1]) != 0:
        raise Exception(f"Expecting first number on first line of '{fileName}' to be 0, but is {parts[1]}")

    ra, dec = map(float, parts[2:4])
    ra *= ANGLE_DEGREE
    dec *= ANGLE_DEGREE
    centerPos = (ra, dec)

    for l in lines:
        if not l.strip() or l.strip().startswith("#"):
            continue

        parts = l.split()
        idNum, ra, dec, a, b, theta, mag, lum = map(float, parts)
        if a != 1 or b != 1 or theta != 0:
            raise Exception(f"Expecting a = 1 (got {a}), b = 1 (got {b}) and theta = 0 (got {theta})")

        ra *= ANGLE_DEGREE
        dec *= ANGLE_DEGREE
        x, y = images.centerOnPosition([ra, dec], centerPos)
        x /= ANGLE_ARCSEC
        y /= ANGLE_ARCSEC
        if not useRADirection:
            x = -x

        velDispFactor = 10**(0.4*(mag0-mag)/vdSlope)
        coreFactor = 10**(0.4*(mag0-mag)/2)
        cutFactor = 10**(0.4*(mag0-mag)*2/slope)
        
        # Better to not mention the hard limits here, are not really used internally but
        # might cause a sanity check to fail after a model refinement
        vdispStr = _getLimitsAndPrior(bsSigma, lambda s: s + " * 1000", cnameStr = f'"cname": "{vdVarName} ; x * {velDispFactor}"', noHard=True)
        cutStr = _getLimitsAndPrior(bsCut, lambda s: s + f" * {arcsecString}", cnameStr = f'"cname": "{cutVarName} ; x * {cutFactor}"', noHard=True)

        if "{" in vdispStr:
            vdispStr += f", # initial range and possibly hard limits are not used, the settings of {vdVarName} will be used, "
        else:
            vdispStr += "(" + vdispStr + f") * {velDispFactor}, # "
        vdispStr += f"factor comes from 10**(0.4*({mag0}-{mag})/{vdSlope})"

        if "{" in cutStr:
            cutStr += ", # initial range and possibly hard limits are not used, the settings of {cutVarName} will be used, "
        else:
            cutStr = "(" + cutStr + f") * {cutFactor}, # "
        cutStr += f"factor comes from 10**(0.4*({mag0} - {mag})*2/{slope})"

        _addBlock(outputLines, {
            "x": f"{x} * {arcsecString},",
            "y": f"{y} * {arcsecString},",
            "angle": "0,",
            "velocitydispersion": vdispStr,
            "coreradius": "(" + coreStr + f") * {coreFactor}, # factor comes from 10**(0.4*({mag0}-{mag})/2)",
            "scaleradius": cutStr,
            "ellipticity": "0",
            "comment": f"Galaxy ID {idNum} from file {fileName}",
        })

    return outputLines

def _getLimitsAndPrior(setting, modifier = lambda s: s, flipMinMax = False, cnameStr = "", noHard=False):
    if int(setting[0]) == 0: # Just a value
        return modifier(setting[1])
    elif int(setting[0]) == 1: # Uniform, use hard limits as specified
        minVal, maxVal = modifier(setting[1]), modifier(setting[2])
        if flipMinMax:
            minVal, maxVal = maxVal, minVal
        if noHard:
            return f'{{ "initmin": {minVal}, "initmax": {maxVal}, {cnameStr} }}'
        return f'{{ "initmin": {minVal}, "initmax": {maxVal}, "hardmin": {minVal}, "hardmax": {maxVal}, {cnameStr} }}'

    elif int(setting[0]) == 3: # Gaussian, use 1 sigma for init range, no hard limits for now
        minVal, maxVal = modifier(f"({setting[1]} - {setting[2]})"), modifier(f"({setting[1]} + {setting[2]})")
        mu, sigma = modifier(setting[1]), modifier(setting[2])
        return f'{{ "initmin": {minVal}, "initmax": {maxVal}, "prior": {{ "type": "gaussian", "params": [ {mu}, {sigma} ] }}, {cnameStr} }}'

def createParametricDescriptionFromLenstoolInput(fileName,
                                                 useRADirection = True,
                                                 lensDistanceString = "Dd",
                                                 arcsecString = "ANGLE_ARCSEC",
                                                 degreeString = "ANGLE_DEGREE",
                                                 kpcString = "DIST_KPC"
                                                 ):
    """This function analyzes a `Lenstool <https://projets.lam.fr/projects/lenstool/wiki>`_
    inversion configuration file (e.g. 'mod.par') and returns a dictionary with
    settings that can be used for a parametric inversion. The :mod:`module <grale.ltparamdesc>`
    itself can also be executed from the command line to generate template inversion files.

    The dictionary that is returned has the following entries:
     - ``"images"``: TODO
     - ``"imagesFileName"``: TODO
     - ``"cosmology"``: TODO
     - ``"positionalUncertainty"``: TODO
     - ``"zd"``: TODO
     - ``"center"``: TODO
     - ``"description"``: TODO

    Arguments:

     - `useRADirection`: TODO
     - `lensDistanceString`: TODO
     - `arcsecString`: TODO
     - `degreeString`: TODO
     - `kpcString`: TODO

    """
    from . import images
    from . import cosmology

    lines = open(fileName, "rt").readlines()

    currentBlock = None
    allBlocks = []

    def _checkUnique(name, optional=False):
        count = 0
        block = None
        blockToDelete = None
        for idx, x in enumerate(allBlocks):
            if x["blockname"] == name:
                count += 1
                block = x
                blockToDelete = idx

        if not optional:
            if count != 1:
                raise Exception(f"Expecting a single entry for '{name}', got {count}")
        else:
            if count > 1:
                raise Exception(f"Expecting at most one entry for '{name}', got {count}")

        if blockToDelete is not None:
            del allBlocks[blockToDelete]
        return block

    lineNum = 0
    for l in lines:
        lineNum += 1

        if not l.strip() or l.strip().startswith("#"):
            continue

        if len(l) == 0:
            raise Exception(f"Unexpected line {lineNum}: " + l.rstrip())

        if l[0] == ' ' or l[0] == '\t':
            if currentBlock is None:
                raise Exception(f"Reading a propery on line {lineNum}, but no block has been set: " + l.rstrip())

            parts = l.strip().split()
            key = parts[0]
            values = parts[1:]
            if key in currentBlock["settings"]:
                raise Exception(f"Setting property {key} multiple times in line {lineNum}:" + l.rstrip())

            if key == "end":
                allBlocks.append(currentBlock)
                currentBlock = None
            else:
                currentBlock["settings"][key] = values
        else:
            if currentBlock is not None:
                raise Exception(f"Already processing a block when reading line {lineNum}: " + l.rstrip())

            parts = l.strip().split()
            blockName = parts[0]
            blockArgs = parts[1:]
            currentBlock = { "blockname": blockName, "blockargs": blockArgs, "settings": { } }

    if not currentBlock or currentBlock["blockname"] != "fini":
        raise Exception("Expecting file to end with 'fini'")

    # Check runmode
    b = _checkUnique("runmode")
    bs = b["settings"]
    if not "inverse" in bs or int(bs["inverse"][0]) != 3:
        raise Exception("No 'inverse' of type 3 in 'runmode'")
    
    if not "reference" in bs or int(bs["reference"][0]) != 3:
        raise Exception("Expecting 'reference' of type 3 in 'runmode'")

    ra, dec = float(bs["reference"][1]), float(bs["reference"][2])
    print(f"Detected center {ra}, {dec}", file=sys.stderr)
    center = [ra*ANGLE_DEGREE, dec*ANGLE_DEGREE]

    b = _checkUnique("image")
    bs = b["settings"]
    if not "forme" in bs or int(bs["forme"][0]) != -1:
        raise Exception("Expecting a 'forme' entry in 'image' with value -1")

    if int(bs["multfile"][0]) != 1:
        raise Exception("Expecting 'multfile' of type 1")

    imageFileName = os.path.join(os.path.dirname(fileName), bs["multfile"][1])
    # TODO: read this here? Or just pass filename?
    imgList = images.readLenstoolInputImagesFile(imageFileName, center, useRADirection)
    print(f"Read {len(imgList)} sources from {imageFileName}", file=sys.stderr)

    positionSigmaForMCMC = float(bs["sigposArcsec"][0])
    print(f"Using positional uncertainty {positionSigmaForMCMC} arcsec in MCMC", file=sys.stderr)
    positionSigmaForMCMC *= ANGLE_ARCSEC

    b = _checkUnique("cosmology")
    bs = b["settings"]

    h = float(bs["H0"][0])/100
    Wm = float(bs["omegaM"][0] if "omegaM" in bs else bs["omega"][0])
    Wv = float(bs["omegaX"][0] if "omegaX" in bs else bs["lambda"][0])
    w = float(bs["wX"][0]) if "wX" in bs else -1
    del bs["H0"]
    if "omegaM" in bs:
        del bs["omegaM"]
    else:
        del bs["omega"]
    if "omegaX" in bs:
        del bs["omegaX"]
    else:
        del bs["lambda"]
    if "wX" in bs:
        del bs["wX"]

    if bs:
        raise Exception(f"Unknown info in cosmological settings: {bs}")
    if abs(Wm + Wv - 1.0) > 1e-8:
        raise Exception("Expecting a flat cosmological model")

    cosm = cosmology.Cosmology(h, Wm, 0, Wv, w)

    # Check for 'grille'
    _checkUnique("grille")
    print("Warning: not using 'grille' settings", file=sys.stderr)

    # Check for 'cline'
    if _checkUnique("cline", True):
        print("Warning: not using 'cline' settings", file=sys.stderr)

    # Check for 'champ'
    _checkUnique("champ")
    print("Warning: not using 'champ' settings", file=sys.stderr)

    # Check for 'cosmolimit'
    if _checkUnique("cosmolimit", True):
        print("Warning: not using 'cosmolimit' settings", file=sys.stderr)

    knownZ = None
    #pprint.pprint(allBlocks, stream=sys.stderr)

    lastBlock = None
    outputLines = """{
    "type": "CompositeLens", "params": [
""".splitlines()

    if useRADirection:
        xString = lambda s: s + f" * (-{arcsecString})"
        xStringComment = " # Flip X to use RA based orientation"
        angString = lambda s: "180 - (" + s + ")"
        angStringComment = " # Use angle supplement for RA based orientation"
    else:
        xString = lambda s: s + f" * {arcsecString}"
        xStringComment = ""
        angString = lambda s: s
        angStringComment = ""

    potfileIndex = 0
    for b in allBlocks:
        bs, ba = b["settings"], b["blockargs"]
        if b["blockname"] == "potentiel":
            if lastBlock:
                _addBlock(outputLines, lastBlock)
                lastBlock = None

            prof = int(bs["profil"][0])
            if prof != 81:
                Exception(f"Can't handle profile {prof}")

            zd = float(bs["z_lens"][0])
            if knownZ is None:
                knownZ = zd
            if knownZ != zd:
                raise Exception(f"Encountered redshift {zd}, but was expecting {knownZ}")

            lastBlock = {
                "x": xString(bs["x_centre"][0]) + "," + xStringComment,
                "y": bs["y_centre"][0] + f" * {arcsecString},",
                "angle": angString(bs["angle_pos"][0]) + "," + angStringComment,
                "coreradius": bs["core_radius"][0] + f" * {arcsecString}," if "core_radius" in bs else bs["core_radius_kpc"][0] + f" * {kpcString}/{lensDistanceString},",
                "scaleradius": bs["cut_radius"][0] + f" * {arcsecString}," if "cut_radius" in bs else bs["cut_radius_kpc"][0] + f" * {kpcString}/{lensDistanceString},",
                "velocitydispersion": bs["v_disp"][0] + " * 1000,",
                "ellipticity": bs["ellipticite"][0] + ",",
                "comment": (" ".join(ba)).strip()
            }

        elif b["blockname"] == "limit":
            if not lastBlock:
                raise Exception("No previous lens info to specify optimization parameters for")

            if "v_disp" in bs:
                lastBlock["velocitydispersion"] = _getLimitsAndPrior(bs["v_disp"], lambda s: s + " * 1000") + ","
            if "angle_pos" in bs:
                lastBlock["angle"] = _getLimitsAndPrior(bs["angle_pos"], angString, useRADirection) + "," + angStringComment
            if "ellipticity" in bs:
                lastBlock["ellipticity"] = _getLimitsAndPrior(bs["ellipticity"]) + ","
            if "x_centre" in bs:
                lastBlock["x"] = _getLimitsAndPrior(bs["x_centre"], xString, useRADirection) + "," + xStringComment
            if "y_centre" in bs:
                lastBlock["y"] = _getLimitsAndPrior(bs["y_centre"], lambda s: s + f" * {arcsecString}") + ","
            if "core_radius" in bs:
                lastBlock["coreradius"] = _getLimitsAndPrior(bs["core_radius"], lambda s: s + f" * {arcsecString}") + ","
            if "cut_radius" in bs:
                lastBlock["scaleradius"] = _getLimitsAndPrior(bs["cut_radius"], lambda s: s + f" * {arcsecString}") + ","

            # Finalize this, don't accept any other modifications
            _addBlock(outputLines, lastBlock)
            lastBlock = None

        elif b["blockname"] == "potfile":

            potfileIndex += 1 # We'll use this to create some variable names
            
            if lastBlock: # Finalize if it exists
                _addBlock(outputLines, lastBlock)
                lastBlock = None

            prof = int(bs["type"][0])
            if prof != 81:
                Exception(f"Can't handle profile type {prof}")

            coreStr = bs["core"][0] + f" * {arcsecString}"
            zd = float(bs["z_lens"][0])
            if knownZ is None:
                knownZ = zd
            if knownZ != zd:
                raise Exception(f"Encountered redshift {zd} in potfile entry, but was expecting {knownZ}")

            mag0 = float(bs["mag0"][0])

            # We'll do this again later to possibly add cname identifiers (not very efficient, I know)
            velDisp0 = _getLimitsAndPrior(bs["sigma"], lambda s: s + f" * 1000")
            cut0 = _getLimitsAndPrior(bs["cut"], lambda s: s + f" * {arcsecString}")
            vd0Name = f"veldispVar{potfileIndex}" if "{" in velDisp0 else None
            cut0Name = f"cut0Var{potfileIndex}" if "{" in cut0 else None
            velDisp0 = _getLimitsAndPrior(bs["sigma"], lambda s: s + f" * 1000", cnameStr="" if not vd0Name else f'"cname": "{vd0Name}"')
            cut0 = _getLimitsAndPrior(bs["cut"], lambda s: s + f" * {arcsecString}", cnameStr="" if not cut0Name else f'"cname": "{cut0Name}"')

            vdSlope = _getLimitsAndPrior(bs["vdslope"])
            if "{" in vdSlope:
                raise Exception("Varying vdslope is currently not supported")
            slope = _getLimitsAndPrior(bs["slope"])
            if "{" in slope:
                raise Exception("Varying slope is currently not supported")

            vdSlope, slope = float(vdSlope), float(slope)

            # Add a component that has no direct effect, but is used to optimize
            # parameters
            _addBlock(outputLines, {
                    "x": "0,",
                    "y": "0,",
                    "angle": "0,",
                    "velocitydispersion": velDisp0 + ",",
                    "coreradius": coreStr + ",",
                    "scaleradius": cut0 + ",",
                    "ellipticity": "0.0000",
                    "comment": "Dummy block to allow potfile parameter changes",
                }, factor=0)

            if int(bs["filein"][0]) != 3:
                raise Exception("Expecting file type 3 for 'filein' in potfile section")

            galFileName = os.path.join(os.path.dirname(fileName), bs["filein"][1])
            outputLines += _processPotfileFile(galFileName, mag0, slope, vdSlope, coreStr, bs["sigma"], bs["cut"], vd0Name, cut0Name,
                                               useRADirection, lensDistanceString, arcsecString, degreeString, kpcString)
        else:
            raise Exception(f"Unknown block type {b['blockname']}")


    if lastBlock:
        _addBlock(outputLines, lastBlock)
        lastBlock = None

    outputLines.append("]}")

    #print("\n".join(outputLines))
    #desc = eval("\n".join(outputLines))
    #paramdesc.analyzeParametricLensDescription(desc, Dd, 0.01)

    return { "images": imgList,
             "imagesFileName": imageFileName,
             "cosmology": cosm,
             "positionalUncertainty": positionSigmaForMCMC,
             "zd": knownZ,
             "center": center,
             "description": "\n".join(outputLines)
           }

def _startFromImgFileName(f, r, useRADirection, outParamFn):
    cosm = r["cosmology"]
    zd = r["zd"]
    Dd = cosm.getAngularDiameterDistance(zd)
    cosmParams = cosm.getParameters()
    centerPos = r["center"]
    imgFn = r["imagesFileName"]
    mcmcUncert = r["positionalUncertainty"]

    print(f"""from grale.all import *
import pickle

zd = {zd}

cosmParams = {cosmParams}
cosm = cosmology.Cosmology(**cosmParams)
cosmology.setDefaultCosmology(cosm)
Dd = cosm.getAngularDiameterDistance(zd)

center = V({centerPos[0]/ANGLE_DEGREE},{centerPos[1]/ANGLE_DEGREE}) * ANGLE_DEGREE
imgList = images.readLenstoolInputImagesFile("{imgFn}", center, {useRADirection})

positionalUncertainty = {mcmcUncert/ANGLE_ARCSEC}*ANGLE_ARCSEC
""", file=f)
    if outParamFn:
        print(f'initialParamDesc = eval(open("{outParamFn}","rt").read()) # Process file as if commands were at this point in the file', file=f)
    else:
        print("initialParamDesc = " + r["description"], file=f)
    print(file=f)

def _startFromLtConfig(f, ltCfgFn, useRADirection):
    print(f"""from grale.all import *
import pickle

r = ltparamdesc.createParametricDescriptionFromLenstoolInput("{ltCfgFn}", useRADirection={useRADirection})

cosm = r["cosmology"]
cosmology.setDefaultCosmology(cosm)

zd = r["zd"]
Dd = cosm.getAngularDiameterDistance(zd)

imgList = r["images"]
positionalUncertainty = r["positionalUncertainty"]

initialParamDesc = eval(r["description"]) # 'description' holds a string, we must interpret it
""", file=f)

def usage():
    print("""
Usage:
    python -m grale.ltparamdesc -in mod.par -out inv.py [-outparam desc.py] [-noRAdir] [-force] [-fromimgfile]
""", file=sys.stderr)
    sys.exit(-1)

def main():
    idx = 1
    ltCfgFn = None
    outFn = None
    outParamFn = None
    useRADir = True
    force = False
    fromImgFileName = False

    try:
        while idx < len(sys.argv):
            opt = sys.argv[idx]
            if opt == "-in":
                ltCfgFn = sys.argv[idx+1]
                idx += 1
            elif opt == "-out":
                outFn = sys.argv[idx+1]
                idx += 1
            elif opt == "-outparam":
                outParamFn = sys.argv[idx+1]
                idx += 1
            elif opt == "-noRAdir":
                useRADir = False
            elif opt == "-force":
                force = True
            elif opt == "-fromimgfile":
                fromImgFileName = True
            else:
                raise Exception(f"Unknown option '{opt}'")

            idx += 1

        if not ltCfgFn:
            raise Exception("Input file needs to be set using '-in'")
        if not outFn:
            raise Exception("Output file needs to be set using '-out'")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        usage()

    try:
        if os.path.exists(outFn) and not force:
            raise Exception(f"Output file '{outFn}' already exists, use '-force' to continue anyway")
        if outParamFn and os.path.exists(outParamFn) and not force:
            raise Exception(f"Output parameters file {outParamFn} already exists, use '-force' to continue anyway")

        r = createParametricDescriptionFromLenstoolInput(ltCfgFn, useRADirection=useRADir)
        if outParamFn:
            open(outParamFn, "wt").write(r["description"])

        with open(outFn, "wt") as f:
            if fromImgFileName:
                _startFromImgFileName(f, r, useRADir, outParamFn)
            else:
                _startFromLtConfig(f, ltCfgFn, useRADir)
                if outParamFn:
                    print(f"Warning: parameter file name '{outParamFn}' will not be used in this inversion mode", file=sys.stderr)

            mcmcGenerations = 5000
            popSize = 512

            print("""# This is a helper function that will be called later. It performs the first
# inversion, not using MCMC but using the genetic algorithm to find a single
# best solution for this parametric description
def initialInversion(imgList, paramDesc):

    # This configures the inversion code to continuously add some random noise
    # to the input positions, to prevent overfitting
    imgListWithUncert = images.addPositionUncertainty(imgList, 0.05*ANGLE_ARCSEC) # TODO: use a different uncertainty here?

    # Create the inversion workspace (the size doesn't really matter for a
    # parametric inversion, but is still a required parameter)
    iws = inversion.InversionWorkSpace(zd, 100*ANGLE_ARCSEC)
    for i in imgListWithUncert:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
        # Since the MCMC version will only use this one, add it so the objectives
        # as well, mainly in case prior terms are used (not used by 'pointimages')
        iws.addImageDataToList(i["imgdata"], i["z"], "bayesstronglensing")

    # Invert!
    startLens, _, _, _, _ = iws.invertParametric(paramDesc, """ + str(popSize) + """, eaType = "JADE",
                                                 fitnessObjectParameters = { "fitness_bayesweaklensing_stronglenssigma": positionalUncertainty }
                                                 )

    startLens.save("startInv.lensdata")

    # This routine will create a modified parametric description, where
    # the initial parameter ranges will be close to the value from the
    # solution that was found. This can then be used in an MCMC run to
    # explore parameter uncertainties
    refinedModel = paramdesc.refineParametricDescription(paramDesc, startLens, 0.001)
    return refinedModel

# This is the second helper function, that will be called with the refined
# parametric description. It uses the Goodman-Weare algorithm (same as 'emcee')
# to explore the parameter space
def mcmcInversion(imgList, paramDesc):

    iws = inversion.InversionWorkSpace(zd, 100*ANGLE_ARCSEC)
    for i in imgList:
        # Here we want to use the logprob based on the deviations in the
        # image plane. This is done by specifying that the images are used
        # in the bayesian code, for the strong lensing contribution.
        iws.addImageDataToList(i["imgdata"], i["z"], "bayesstronglensing")

    # Well run the algorithm for several generation, half of which will be
    # ignored as burn-in samples
    totalGenerations = """ + str(mcmcGenerations) + """
    burnInGenerations = totalGenerations//2

    # We'll also specify the file to which the samples are written
    mcmcParams = { "burningenerations": burnInGenerations, "samplesfilename": "samples.dat" }

    # This specifies that all strong lensing images are considered to have an
    # uncertainty
    fitParams = { "fitness_bayesweaklensing_stronglenssigma": positionalUncertainty }

    # Run the MCMC algorithm! It will move several 'walkers' around, for several generations.
    lens, _, _, reducedParamInfo, _ = iws.invertParametric(paramDesc, """ + str(popSize) + """, maximumGenerations = totalGenerations, eaType = "MCMC", 
                                                 geneticAlgorithmParameters=mcmcParams, fitnessObjectParameters=fitParams,
                                                 clampToHardLimits=True)

    lens.save("mcmcInv.lensdata")

    # To interpret the samples that are written to file, we also write out the
    # information about the parameters that are being varied. These parameters
    # typically use some find of scale factor, to allow for 32-bit floating point
    # calculations, and it's those scaled values that are written to the file.
    # The 'paramInfo' list contains the scale factors to convert the values from
    # file to the real values.
    pickle.dump(reducedParamInfo, open("paraminfo.dat", "wb"))

# Run the initial inversion to get a good starting point for the MCMC part
newParamDesc = initialInversion(imgList, initialParamDesc)
# Do the actual MCMC inversion
mcmcInversion(imgList, newParamDesc)

""", file=f)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(-1)

if __name__ == "__main__":
    main()

