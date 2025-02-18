"""This is a module with a helper function to import parametric inversion
settings from a Lenstool configuration file. The module itself can also
be executed on the command line to generate template inversion scripts
based on a Lenstool configuration.

Usage::

    python -m grale.ltparamdesc -in mod.par (-out inv.py | -exec)\\
          [-outparam desc.py] [-noRAdir] [-force] [-fromimgfile] \\
          [-popsize 512] [-mcmcgen 5000] [-outfnsuffix suff] \\
          [-refinefactor 0.0001] [-debug] [-nostricimages]

Arguments:

 - `-in mod.par`: this refers to the Lenstool input file that should be
   processed.
 - `-out inv.py`: this specifies the name of the main Python script that
   is to be generated. The script will contain code to setup the inversion
   and run it in two steps: in the first part, a good starting point for
   all the parameters will be determined, using a GA-like approach. In
   the second step the parameters around this point will be explored using
   an MCMC technique to reveal information about the uncertainties of the
   parameters.
 - `-exec`: if set, the `-out` option may be omitted and the resulting
   script will be executed
 - `-outparam desc.py`: if specified, this optional parameter gives the
   name of the file that will contain the parametric description that
   has been derived from the Lenstool input file. See the :mod:`paramdesc <grale.paramdesc>`
   module for an explanation of such a description.
 - `-noRAdir`: by default, the coordinates in the input file will be
   converted so that they follow the orientation of the RA axis (so this axis
   points left). To mirror the coordinates you can specify this flag.
 - `-force`: by default, the program will halt instead of overwriting one
   of the output files. To disable this safety you can set this parameter.
 - `-fromimgfile`: the script that's generated will again read the same
   file as specified using the `-in` option. You can also make it so that
   the script only reads the images file, by setting this command line
   option. If no `-outparam` is specified, then the entire parametric
   description is stored in the output script. If the `-outparam` is
   specified, the generated script is smaller and will read the model
   from that file.
 - `-popsize 512`: sets the population size for initial inversion and MCMC
   parts in the generated script. Defaults to 512
 - `-mcmcgen 5000`: sets the number of generations/steps to perform during
   the MCMC exploration, half of this will be considered as burn-in. Defaults
   to 5000.
 - `-outfnsuffix suff`: for the output file names in the script, use this
   as the suffix (default is empty).
 - `-refinefactor 0.0001`: when the initial inversion has been found, the
   MCMC phase will start from a narrow region around that model's parameters.
   This specifies the fraction of the parameter change that's used for the
   initialization of a new set of trial solutions.
 - `-debug`: if present, more debug output will be shown.
 - `-nostrictimages`: if set, the input images file is not required to
   start with '#REFERENCE 0'
"""
from .constants import *
from . import cosmology
from . import images
from . import inversion
from . import paramdesc
import sys
import os
from io import StringIO
import traceback
import pickle

# "x", "y", "angle", "velocitydispersion", "coreradius", "scaleradius", "ellipticity" entries
def _addBlock(outputLines, block, factor=1):

    lensName = "LTPIEMDLens"
    try:
        ellStr = block["ellipticity"].strip()
        if ellStr[-1] == ",":
            ellStr = ellStr[:-1]
        ellipt = float(ellStr)
        if ellipt == 0:
            lensName = "LTPIMDLens"
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
            "type": "{lensName}",
            "params": {{
                "velocitydispersion": {block['velocitydispersion']}
                "coreradius": {block['coreradius']}
                "scaleradius": {block['scaleradius']}
"""
    if lensName == "LTPIEMDLens":
        o += f"""                "ellipticity": {block['ellipticity']}
"""
    o += """            }
        }
    },"""
    for l in o.splitlines():
        outputLines.append(l)

def _processPotfileFile(fileName, mag0, slope, vdSlope, coreStr, bsSigma, bsCut, isCutArcsec, vdVarName, cutVarName,
                        useRADirection, lensDistanceString, arcsecString, degreeString, kpcString,
                        useRelativeRaDec):

    outputLines = []

    lines = open(fileName, "rt").readlines()

    if not useRelativeRaDec:
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
        #if a != 1 or b != 1 or theta != 0:
            #raise Exception(f"Expecting a = 1 (got {a}), b = 1 (got {b}) and theta = 0 (got {theta})")
            #print(f"Warning in potfile {fileName}: Expecting a = 1 (got {a}), b = 1 (got {b}) and theta = 0 (got {theta})")

        ellipticity = (a**2-b**2)/(a**2+b**2)

        if useRelativeRaDec:
            x = ra*ANGLE_ARCSEC
            y = dec*ANGLE_ARCSEC
            if useRADirection: # TODO: is this correct?
                x = -x
        else:
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
        vdispStr = _getLimitsAndPrior(bsSigma, lambda s: s + " * 1000", cnameStr = f'"cname": "{vdVarName} ; x * {velDispFactor}"', noHard=True, noPrior=True)

        cutLambda = (lambda s: s + f" * {arcsecString}") if isCutArcsec else (lambda s: s + f" * {kpcString}/{lensDistanceString}")
        cutStr = _getLimitsAndPrior(bsCut, cutLambda, cnameStr = f'"cname": "{cutVarName} ; x * {cutFactor}"', noHard=True, noPrior=True)

        if "{" in vdispStr:
            vdispStr += f", # initial range and possibly hard limits are not used, the settings of {vdVarName} will be used, "
        else:
            vdispStr = "(" + vdispStr + f") * {velDispFactor}, # "
        vdispStr += f"factor comes from 10**(0.4*({mag0}-{mag})/{vdSlope})"

        if "{" in cutStr:
            cutStr += ", # initial range and possibly hard limits are not used, the settings of {cutVarName} will be used, "
        else:
            cutStr = "(" + cutStr + f") * {cutFactor}, # "
        cutStr += f"factor comes from 10**(0.4*({mag0} - {mag})*2/{slope})"

        _addBlock(outputLines, {
            "x": f"{x} * {arcsecString},",
            "y": f"{y} * {arcsecString},",
            "angle": f"{theta}," if not useRADirection else f"(180-{theta}), # Use angle supplement for RA based orientation",
            "velocitydispersion": vdispStr,
            "coreradius": "(" + coreStr + f") * {coreFactor}, # factor comes from 10**(0.4*({mag0}-{mag})/2)",
            "scaleradius": cutStr,
            "ellipticity": f"{ellipticity}",
            "comment": f"Galaxy ID {idNum} from file {fileName}",
        })

    return outputLines

def _getLimitsAndPrior(setting, modifier = lambda s: s, flipMinMax = False, cnameStr = "", noHard=False, noPrior=False):
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

        if noPrior:
            # This is set if we're using a cname to derive this variable's value from another one
            # In that case, only the original variable should have prior information set
            return f'{{ "initmin": {minVal}, "initmax": {maxVal}, {cnameStr} }}'

        return f'{{ "initmin": {minVal}, "initmax": {maxVal}, "prior": {{ "type": "gaussian", "params": [ {mu}, {sigma} ] }}, {cnameStr} }}'

    raise Exception(f"Unrecognized block setting '{setting}'")

def createParametricDescriptionFromLenstoolInput(fileName,
                                                 useRADirection,
                                                 parseImages = True,
                                                 strictImages = True,
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

     - ``"images"``: this is a list of the images, basically what the
       :func:`readLenstoolInputImagesFile <grale.images.readLenstoolInputImagesFile>`
       returns for the image file specified in the Lenstool script. Is ``None``
       in case the `parseImages` flag was set to ``False``.
     - ``"imagesFileName"``: the name of the images file that is stored in the Lenstool
       script.
     - ``"cosmology"``: an instance of :class:`Cosmology <grale.cosmology.Cosmology>`,
       based on the settings that were read from the input file.
     - ``"positionalUncertainty"``: the uncertainly that is applicable to the input images
     - ``"zd"``: the redshift of the lens
     - ``"center"``: this should be used as the center of the coordinate system, can
       be used in a call to :func:`readLenstoolInputImagesFile <grale.images.readLenstoolInputImagesFile>`
       for example.
     - ``"description"``: a string containing the parametric description (see module 
       :mod:`paramdesc <grale.paramdesc>`). You'll probably need to '`eval <https://docs.python.org/3/library/functions.html#eval>`_'
       this to be able to actually use it, but the string representation is much more
       readable.

    Arguments:

     - `fileName`: name of the Lenstool input file that needs to be analyzed
     - `useRADirection`: if ``True``, coordinates in the Lenstool file are converted
       so that they use the RA direction (axis points left). To use mirrored coordinates
       (as Lenstool itself does) you can set this to ``False``
     - `parseImages`: can be set to ``False`` if the specified images file is not
       processed further, only the file name will be stored.
     - `strictImages`: by default, the input images file needs to start with
       '#REFERENCE 0', by setting this flag to ``False`` this requirement is relaxed.
     - `lensDistanceString`: sometimes a conversion may be needed that involves the
       lens distance. This string is used in the generated output for this.
     - `arcsecString`: string to be used when expressing something in arcsec
     - `degreeString`: string to be used when expressing something in degrees
     - `kpcString`: string to be used when expressing something in kpc.

    """

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
            if blockName == "cosmologie":
                blockName = "cosmology"
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
    if parseImages:
        imgList = images.readLenstoolInputImagesFile(imageFileName, center, useRADirection, strictImages)
        print(f"Read {len(imgList)} sources from {imageFileName}", file=sys.stderr)
    else:
        imgList = None
        print(f"Skipping reading images file {imageFileName}")

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

            if "ellipticite" in bs and "ellipticity" in bs:
                raise Exception("Both 'ellipticity' and 'ellipticite' present in limit block")
        
            # It seems that if 0 is the setting in the limit block, then the value
            # from the previous 'potentiel' block should be used
            def _checkPresentAndVariable(n):
                if n in bs and int(bs[n][0]) != 0:
                    return True
                return False

            if _checkPresentAndVariable("v_disp"):
                lastBlock["velocitydispersion"] = _getLimitsAndPrior(bs["v_disp"], lambda s: s + " * 1000") + ","
            if _checkPresentAndVariable("angle_pos"):
                lastBlock["angle"] = _getLimitsAndPrior(bs["angle_pos"], angString, useRADirection) + "," + angStringComment
            if _checkPresentAndVariable("ellipticity"):
                lastBlock["ellipticity"] = _getLimitsAndPrior(bs["ellipticity"]) + ","
            if _checkPresentAndVariable("ellipticite"):
                lastBlock["ellipticity"] = _getLimitsAndPrior(bs["ellipticite"]) + ","
            if _checkPresentAndVariable("x_centre"):
                lastBlock["x"] = _getLimitsAndPrior(bs["x_centre"], xString, useRADirection) + "," + xStringComment
            if _checkPresentAndVariable("y_centre"):
                lastBlock["y"] = _getLimitsAndPrior(bs["y_centre"], lambda s: s + f" * {arcsecString}") + ","
            if _checkPresentAndVariable("core_radius"):
                lastBlock["coreradius"] = _getLimitsAndPrior(bs["core_radius"], lambda s: s + f" * {arcsecString}") + ","
            if _checkPresentAndVariable("core_radius_kpc"):
                lastBlock["coreradius"] = _getLimitsAndPrior(bs["core_radius_kpc"], lambda s: s + f" * {kpcString}/{lensDistanceString}") + ","
            if _checkPresentAndVariable("cut_radius"):
                lastBlock["scaleradius"] = _getLimitsAndPrior(bs["cut_radius"], lambda s: s + f" * {arcsecString}") + ","
            if _checkPresentAndVariable("cut_radius_kpc"):
                lastBlock["scaleradius"] = _getLimitsAndPrior(bs["cut_radius_kpc"], lambda s: s + f" * {kpcString}/{lensDistanceString}") + ","

            # Finalize this, don't accept any other modifications
            _addBlock(outputLines, lastBlock)
            lastBlock = None

        elif b["blockname"] == "potfile":

            potfileIndex += 1 # We'll use this to create some variable names
            
            if lastBlock: # Finalize if it exists
                _addBlock(outputLines, lastBlock)
                lastBlock = None

            prof = int(bs["type"][0]) if "type" in bs else 81 # Defaults to 81 I think
            if prof != 81:
                Exception(f"Can't handle profile type {prof}")

            coreStr = (bs["core"][0] + f" * {arcsecString}") if "core" in bs else (bs["corekpc"][0] + f" * {kpcString}/{lensDistanceString}")

            if "z_lens" in bs:
                zd = float(bs["z_lens"][0])
                if knownZ is None:
                    knownZ = zd
                if knownZ != zd:
                    raise Exception(f"Encountered redshift {zd} in potfile entry, but was expecting {knownZ}")

            mag0 = float(bs["mag0"][0])

            # We'll do this again later to possibly add cname identifiers (not very efficient, I know)
            velDisp0 = _getLimitsAndPrior(bs["sigma"], lambda s: s + f" * 1000")
            cut0 = _getLimitsAndPrior(bs["cut"], lambda s: s + f" * {arcsecString}") if "cut" in bs else _getLimitsAndPrior(bs["cutkpc"], lambda s: s + f" * {kpcString}/{lensDistanceString}")
            vd0Name = f"veldispVar{potfileIndex}" if "{" in velDisp0 else None
            cut0Name = f"cut0Var{potfileIndex}" if "{" in cut0 else None
            velDisp0 = _getLimitsAndPrior(bs["sigma"], lambda s: s + f" * 1000", cnameStr="" if not vd0Name else f'"cname": "{vd0Name}"')
            
            cutCname = "" if not cut0Name else f'"cname": "{cut0Name}"' 
            cut0 = _getLimitsAndPrior(bs["cut"], lambda s: s + f" * {arcsecString}", cnameStr=cutCname) if "cut" in bs else _getLimitsAndPrior(bs["cutkpc"], lambda s: s + f" * {kpcString}/{lensDistanceString}", cnameStr=cutCname)

            vdSlope = _getLimitsAndPrior(bs["vdslope"] if "vdslope" in bs else [0, "4.0"])
            if "{" in vdSlope:
                raise Exception("Varying vdslope is currently not supported")
            slope = _getLimitsAndPrior(bs["slope"] if "slope" in bs else [0, "4.0"])
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


            potFileIn = int(bs["filein"][0])
            if not potFileIn in [1, 3]:
                raise Exception("Expecting file type 1 or 3 for 'filein' in potfile section")
            useRelativeRaDec = True if potFileIn == 1 else False

            galFileName = os.path.join(os.path.dirname(fileName), bs["filein"][1])

            outputLines += _processPotfileFile(galFileName, mag0, slope, vdSlope, coreStr, bs["sigma"], bs["cut"] if "cut" in bs else bs["cutkpc"],
                                               "cut" in bs, vd0Name, cut0Name,
                                               useRADirection, lensDistanceString, arcsecString, degreeString, kpcString,
                                               useRelativeRaDec)
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

def _startFromImgFileName(f, r, useRADirection, outParamFn, debugCode, strictImages):
    cosm = r["cosmology"]
    zd = r["zd"]
    Dd = cosm.getAngularDiameterDistance(zd)
    cosmParams = cosm.getParameters()
    centerPos = r["center"]
    imgFn = r["imagesFileName"]
    mcmcUncert = r["positionalUncertainty"]

    print(f"""from grale.all import *
import pickle
{debugCode}
zd = {zd}

cosmParams = {cosmParams}
cosm = cosmology.Cosmology(**cosmParams)
cosmology.setDefaultCosmology(cosm)
Dd = cosm.getAngularDiameterDistance(zd)

center = V({centerPos[0]/ANGLE_DEGREE},{centerPos[1]/ANGLE_DEGREE}) * ANGLE_DEGREE
imgList = images.readLenstoolInputImagesFile("{imgFn}", center, {useRADirection}, strict={strictImages})

positionalUncertainty = {mcmcUncert/ANGLE_ARCSEC}*ANGLE_ARCSEC
""", file=f)
    if outParamFn:
        print(f'initialParamDesc = eval(open("{outParamFn}","rt").read()) # Process file as if commands were at this point in the file', file=f)
    else:
        print("initialParamDesc = " + r["description"], file=f)
    print(file=f)

def _startFromLtConfig(f, ltCfgFn, useRADirection, debugCode, strictImages):
    print(f"""from grale.all import *
import pickle
{debugCode}
r = ltparamdesc.createParametricDescriptionFromLenstoolInput("{ltCfgFn}", useRADirection={useRADirection}, strictImages={strictImages})

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
    python -m grale.ltparamdesc -in mod.par (-out inv.py | -exec) [-outparam desc.py] [-noRAdir] [-force] [-fromimgfile] \\
                                [-mcmcgen 5000] [-popsize 512] [-outfnsuffix] [-refinefactor 0.0001] [-debug] \\
                                [-nostrictimages]
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
    execScript = False
    mcmcGenerations = 5000
    popSize = 512
    outFnSuffix = ""
    refineFactor = 0.0001
    strictImages = True
    debug = False

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
            elif opt == "-exec":
                execScript = True
            elif opt == "-mcmcgen":
                mcmcGenerations = int(sys.argv[idx+1])
                idx += 1
            elif opt == "-popsize":
                popSize = int(sys.argv[idx+1])
                idx += 1
            elif opt == "-outfnsuffix":
                outFnSuffix = sys.argv[idx+1]
                idx += 1
            elif opt == "-refinefactor":
                refineFactor = float(sys.argv[idx+1])
                idx += 1
            elif opt == "-debug":
                debug = True
            elif opt == "-nostrictimages":
                strictImages = False
            else:
                raise Exception(f"Unknown option '{opt}'")

            idx += 1

        if not ltCfgFn:
            raise Exception("Input file needs to be set using '-in'")
        if not execScript and not outFn:
            raise Exception("Output file needs to be set using '-out' if '-exec' is not set")
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        usage()

    debugCode = ""
    if debug:
        debugCode = """
# May be helpful for debugging
inverters.debugOutput = True
inverters.debugDirectStderr = True
"""

    try:
        if outFn and os.path.exists(outFn) and not force:
            raise Exception(f"Output file '{outFn}' already exists, use '-force' to continue anyway")
        if outParamFn and os.path.exists(outParamFn) and not force:
            raise Exception(f"Output parameters file {outParamFn} already exists, use '-force' to continue anyway")

        r = createParametricDescriptionFromLenstoolInput(ltCfgFn, useRADirection=useRADir, parseImages=(not fromImgFileName), strictImages=strictImages)
        if outParamFn:
            open(outParamFn, "wt").write(r["description"])

        stringOutput = StringIO()
        if fromImgFileName:
            _startFromImgFileName(stringOutput, r, useRADir, outParamFn, debugCode, strictImages)
        else:
            _startFromLtConfig(stringOutput, ltCfgFn, useRADir, debugCode, strictImages)
            if outParamFn:
                print(f"Warning: parameter file name '{outParamFn}' will not be used in this inversion mode", file=sys.stderr)

        print(f"""# This is a helper function that will be called later. It performs the first
# inversion, not using MCMC but using the genetic algorithm to find a single
# best solution for this parametric description
def initialInversion(zd, imgList, paramDesc, positionalUncertainty):

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
    startLens, _, _, _, _ = iws.invertParametric(paramDesc, {popSize}, eaType = "JADE",
                                                 fitnessObjectParameters = {{ "fitness_bayesweaklensing_stronglenssigma": positionalUncertainty }},
                                                 useImagePositionRandomization=True, # Needs to be set for the addPositionUncertainty to have effect
                                                 )

    startLens.save("startInv{outFnSuffix}.lensdata")

    # This routine will create a modified parametric description, where
    # the initial parameter ranges will be close to the value from the
    # solution that was found. This can then be used in an MCMC run to
    # explore parameter uncertainties
    refinedModel = paramdesc.refineParametricDescription(paramDesc, startLens, {refineFactor})
    return refinedModel

# This is the second helper function, that will be called with the refined
# parametric description. It uses the Goodman-Weare algorithm (same as 'emcee')
# to explore the parameter space
def mcmcInversion(zd, imgList, paramDesc, positionalUncertainty):

    iws = inversion.InversionWorkSpace(zd, 100*ANGLE_ARCSEC)
    for i in imgList:
        # Here we want to use the logprob based on the deviations in the
        # image plane. This is done by specifying that the images are used
        # in the bayesian code, for the strong lensing contribution.
        iws.addImageDataToList(i["imgdata"], i["z"], "bayesstronglensing")

    # Well run the algorithm for several generation, half of which will be
    # ignored as burn-in samples
    totalGenerations = {mcmcGenerations}
    burnInGenerations = totalGenerations//2

    # We'll also specify the file to which the samples are written
    mcmcParams = {{ "burningenerations": burnInGenerations, "samplesfilename": "samples{outFnSuffix}.dat" }}

    # This specifies that all strong lensing images are considered to have an
    # uncertainty
    fitParams = {{ "fitness_bayesweaklensing_stronglenssigma": positionalUncertainty }}

    # Run the MCMC algorithm! It will move several 'walkers' around, for several generations.
    lens, _, _, reducedParamInfo, _ = iws.invertParametric(paramDesc, {popSize}, maximumGenerations = totalGenerations, eaType = "MCMC", 
                                                 geneticAlgorithmParameters=mcmcParams, fitnessObjectParameters=fitParams,
                                                 clampToHardLimits=True)

    lens.save("mcmcInv{outFnSuffix}.lensdata")

    # To interpret the samples that are written to file, we also write out the
    # information about the parameters that are being varied. These parameters
    # typically use some find of scale factor, to allow for 32-bit floating point
    # calculations, and it's those scaled values that are written to the file.
    # The 'paramInfo' list contains the scale factors to convert the values from
    # file to the real values.
    pickle.dump(reducedParamInfo, open("paraminfo{outFnSuffix}.dat", "wb"))

# Run the initial inversion to get a good starting point for the MCMC part
newParamDesc = initialInversion(zd, imgList, initialParamDesc, positionalUncertainty)
# Do the actual MCMC inversion
mcmcInversion(zd, imgList, newParamDesc, positionalUncertainty)

""", file=stringOutput)

        finalScript = stringOutput.getvalue()
        if outFn:
            with open(outFn, "wt") as f:
                f.write(finalScript)

        if execScript:
            print("Executing final script:")
            exec(finalScript)

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(-1)

if __name__ == "__main__":
    main()

