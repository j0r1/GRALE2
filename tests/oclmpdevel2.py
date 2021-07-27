from grale.all import *
import pprint
import pyopencl as cl

# Dmatrix[i,j]:
# i = 0..numPlanes (inclusive)
# j = 0..(numPlanes-1)
# i = 0 = Dist(z=0, zd_(j+1))
# i != 0 = Dist(zd_i,zd_(j+1))
# Dsrc[0] = Dist(0, z_s)
# Dsrc[i] = Dist(zd_i, z_s)

def checkPlanes(z_ds):
    numPlanes = len(z_ds)
    for j in range(1,numPlanes):
        if z_ds[j] <= z_ds[j-1]:
            raise Exception("Lens plane redshifts must be increasing")
    return numPlanes

def getDistanceMatrix(z_ds, cosm):

    numPlanes = checkPlanes(z_ds)
    Dmatrix = np.ones(shape=(numPlanes+1, numPlanes), dtype=np.double) * float("NaN")

    for j in range(numPlanes):
        Dmatrix[0,j] = cosm.getAngularDiameterDistance(z_ds[j])/DIST_MPC

    # There are some double calculations now
    for i in range(1, numPlanes+1):
        for j in range(i, numPlanes):
            Dmatrix[i,j] = cosm.getAngularDiameterDistance(z_ds[i-1], z_ds[j])/DIST_MPC

    return Dmatrix.astype(np.float32)
    
def getDistanceVector(z_s, z_ds, cosm):
    numPlanes = checkPlanes(z_ds)
    Dsrc = np.empty(shape=(numPlanes+1,), dtype=np.double)
    Dsrc[0] = cosm.getAngularDiameterDistance(z_s)/DIST_MPC
    usedPlanes = 0
    for i in range(0, numPlanes):
        if z_s > z_ds[i]:
            usedPlanes += 1
            Dsrc[i+1] = cosm.getAngularDiameterDistance(z_ds[i], z_s)/DIST_MPC
        else:
            Dsrc[i+1] = float("NaN")

    return Dsrc, usedPlanes

# TODO: same types should be consecutive, will be much more performant
def getAlphaCodeForPlane(functionName, plane, neededBasisFunctionCode):

    code = """
float2 {functionName}(float2 theta, __global const int *pIntParams, __global const float *pFloatParams, __global const float *pWeights,
                      __global const float *pCenters, float scalableFunctionWeight)
{{
    float2 alpha = (float2)(0, 0);
    //printf("{functionName}\\n");
""".format(functionName=functionName)

    def addCodeForLens(functionName, num, iCnt, fCnt, usePlaneScale):
        code = ""

        if num != 1:
            code += """
    for (int i = 0 ; i < {Nlenses} ; i++)
    """.format(Nlenses = num)

        code += """
    {{
        const float2 center = (float2)(pCenters[0], pCenters[1]);
        //printf("center = %g %g\\n", center.x, center.y);
        LensQuantities l = {lensFunctionName}(theta-center, pIntParams, pFloatParams);
        const float w = *pWeights;
 
        pWeights++;
        pCenters += 2;
""".format(lensFunctionName=functionName)
        if usePlaneScale:
            code += """
        alpha.x += w*l.alphaX*scalableFunctionWeight;
        alpha.y += w*l.alphaY*scalableFunctionWeight;
"""
        else:
            code += """
        alpha.x += w*l.alphaX;
        alpha.y += w*l.alphaY;
"""

        if iCnt:
            code += "        pIntParams += {};\n".format(iCnt)
        if fCnt:
            code += "        pFloatParams += {};\n".format(fCnt)
        code += """
    }
"""
        return code

    def checkFunction(t, fnName, fnCode):
        # For a compositelens, we'll add the code later (may need other sublenses or more recursion)
        if t == lenses.CompositeLens:
            fnCode = ""

        if not t in neededBasisFunctionCode:
            neededBasisFunctionCode[t] = {
                "functionname": fnName,
                "functioncode": fnCode,
            }

    prevBfType, prevBfCounts, prevIntCount, prevFloatCount = None, 0, None, None
    for bf,center in plane["scaledlenses"]:
        t = type(bf)
        fnName, fnCode = bf.getCLProgram(False, False)
        checkFunction(t, fnName, fnCode)

        iCnt, fCnt = bf.getCLParameterCounts()
        #print("Type", t, iCnt, fCnt)

        if prevBfType != t or (prevBfType == t and (iCnt != prevIntCount or fCnt != prevFloatCount )):
            #print(prevBfType, t, iCnt, prevIntCount, fCnt, prevFloatCount)
            if prevBfCounts:
                code += addCodeForLens(neededBasisFunctionCode[prevBfType]["functionname"], prevBfCounts, prevIntCount, prevFloatCount, True)
            prevBfType, prevBfCounts, prevIntCount, prevFloatCount = t, 1, iCnt, fCnt
        else:
            prevBfCounts += 1

    if prevBfCounts:
        code += addCodeForLens(neededBasisFunctionCode[prevBfType]["functionname"], prevBfCounts, prevIntCount, prevFloatCount, True)

    if plane["unscaledlens"]:
        bf, center = plane["unscaledlens"]
        fnName, fnCode = bf.getCLProgram(False, False)
        checkFunction(type(bf), fnName, fnCode)

        iCnt, fCnt = bf.getCLParameterCounts()
        code += addCodeForLens(fnName, 1, iCnt, fCnt, False)

    code += """
    return alpha;
}""";
    return code

def checkCompositeLenses(lensPlanes):

    totalMaxRecurs = 0
    totalSubCode = { }
    totalSubLensArray = [ ]

    for plane in lensPlanes:
        planeLenses = [ bf for bf, center in plane["scaledlenses"] ]
        if plane["unscaledlens"]:
            planeLenses.append(plane["unscaledlens"][0])

        for bf in planeLenses:
            if type(bf) != lenses.CompositeLens:
                continue

            maxRecurs, subCode, subLensArray = bf.findCLSubroutines(False, False)

            totalMaxRecurs = max(totalMaxRecurs, maxRecurs)
            for x in subCode:
                totalSubCode[x] = subCode[x]

            if not totalSubLensArray:
                for i in subLensArray:
                    totalSubLensArray.append(i)
            else:
                assert(len(subLensArray) == len(totalSubLensArray))
                for i in range(len(subLensArray)):
                    if not totalSubLensArray[i]:
                        totalSubLensArray[i] = subLensArray[i]

    return totalMaxRecurs, totalSubCode, totalSubLensArray

def getMultiPlaneTraceCode(lensPlanes):

    compRecurs, compSubroutineCode, compSubLensArray = checkCompositeLenses(lensPlanes)

    maxPlanes = len(lensPlanes)

    alphaCode = ""
    neededBasisFunctionCode = { }

    code = """

#define MAXPLANES {maxPlanes}

float2 getAlpha(int lpIdx, float2 theta, __global const int *pAllIntParams, __global const float *pAllFloatParams,
                __global const float *pAllWeights, __global const float *pAllCenters, __global const int *pPlaneIntParamOffsets,
                __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets,
                __global const float *pScalableFunctionWeights)
{{

""".format(maxPlanes=maxPlanes)
    
    for i in range(maxPlanes):
        alphaCode += getAlphaCodeForPlane("getAlpha_{}".format(i), lensPlanes[i], neededBasisFunctionCode)
        code += """
    if (lpIdx == {i})
        return getAlpha_{i}(theta, pAllIntParams + pPlaneIntParamOffsets[{i}], pAllFloatParams + pPlaneFloatParamOffsets[{i}],
                            pAllWeights + pPlaneWeightOffsets[{i}], pAllCenters + pPlaneWeightOffsets[{i}]*2, pScalableFunctionWeights[{i}]);
""".format(i = i);

    code += """

    return (float2)(0.0f/0.0f, 0.0f/0.0f);
}

float2 multiPlaneTrace(float2 theta, int numPlanes, __global const float *Dsrc, __global const float *Dmatrix,
                __global const int *pAllIntParams, __global const float *pAllFloatParams, __global const float *pAllWeights,
                __global const float *pAllCenters, __global const int *pPlaneIntParamOffsets,
                __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets,
                __global const float *pScalableFunctionWeights)
{
    float2 T[MAXPLANES];
    float2 alphas[MAXPLANES];

    for (int j = 0 ; j < numPlanes ; j++)
    {
        T[j] = theta;
        for (int i = 0 ; i < j-1 ; i++)
            T[j] -= (Dmatrix[(i+1)*MAXPLANES + j]/Dmatrix[0 + j]) * alphas[i];
        const int i = j-1;
        if (i >= 0)
        {
            const float2 alpha = getAlpha(i, T[i], pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters, pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets, pScalableFunctionWeights);
            alphas[i] = alpha;
	        const float dfrac = Dmatrix[(i+1)*MAXPLANES + j]/Dmatrix[0 + j];

            T[j] -= dfrac * alpha;
        }
    }

    float2 beta = theta;
    for (int i = 0 ; i < numPlanes-1 ; i++)
        beta -= (Dsrc[i+1]/Dsrc[0]) * alphas[i];

    const int i = numPlanes-1;
    if (i >= 0)
    {
        const float2 alpha = getAlpha(i, T[i], pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters, pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets, pScalableFunctionWeights);
        const float dfrac = Dsrc[i+1]/Dsrc[0];
        beta -= dfrac * alpha;
    }
    return beta;
}
"""

    # Merge composite lens subroutines with other subroutines
    allSubRoutines = compSubroutineCode.copy()
    for t in neededBasisFunctionCode:
        n, c = neededBasisFunctionCode[t]["functionname"], neededBasisFunctionCode[t]["functioncode"]
        allSubRoutines[n] = c

    # Add regular subroutines
    allCode = ""
    for n in allSubRoutines:
        allCode += allSubRoutines[n]

    # Add the CompositeLens code if needed
    if compSubLensArray:
        n, c = lenses.CompositeLens.getCompositeCLProgram(compSubLensArray, compRecurs, False, False)
        allCode += c

    # Add the code for the different lens planes
    allCode += alphaCode

    # Add the overall trace code
    allCode += code

    return allCode

def getMultiPlaneOCLProgram(lensplanes):

    code = lensplanes[0]["scaledlenses"][0][0].getCLLensQuantitiesStructure(False, False) + getMultiPlaneTraceCode(lensplanes) + """
__kernel void calculateBetas(const int numPoints, __global const float *pThetas, __global float *pBetas, 
                                __global const int *pNumPlanes, __global const float *DsrcAll, __global const float *Dmatrix,
                                __global const int *pAllIntParams, __global const float *pAllFloatParams,
                                __global const float *pAllWeights, __global const float *pAllCenters,
                                __global const int *pPlaneIntParamOffsets,
                                __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets,
                                __global const float *pScalableFunctionWeights)
{
    const int i = get_global_id(0);
    if (i >= numPoints)
        return;

    const float2 theta = (float2)(pThetas[i*2+0], pThetas[i*2+1]);
    // Each Dsrc is vector of MAXPLANES+1 length
    __global const float *Dsrc = DsrcAll + (MAXPLANES+1)*i;
    const float2 beta = multiPlaneTrace(theta, pNumPlanes[i], Dsrc, Dmatrix,
                                        pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters,
                                        pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets,
                                        pScalableFunctionWeights);

    pBetas[i*2+0] = beta.x;
    pBetas[i*2+1] = beta.y;
}
"""
    return code

def getOpenCLData(cosm, zss, zds, lensplanes, angularScale):

    numPoints = len(zss)
    assert(numPoints == len(zss))

    maxPlanes = len(lensplanes)
    DsrcAll = np.ones((numPoints, maxPlanes+1), dtype=np.float32)*float("NaN")
    allNumPlanes = []
    for i in range(len(zss)):
        Dsrc, usedPlanes = getDistanceVector(zss[i], zds, cosm)
        DsrcAll[i,:len(Dsrc)] = np.array(Dsrc)
        allNumPlanes.append(usedPlanes)

    allNumPlanes = np.array(allNumPlanes, dtype=np.int32)
    Dmatrix = getDistanceMatrix(zds, cosm)

    intParams = []
    floatParams = []
    weights = []
    centers = []
    intParamOffsets = [ ]
    floatParamOffsets = [ ]
    weightOffsets = [ ]
    planeWeights = [ ]

    for plane in lensplanes:
        intParamOffsets.append(len(intParams))
        floatParamOffsets.append(len(floatParams))
        weightOffsets.append(len(weights))
        planeWeights.append(1.0)

        allLenses = plane["scaledlenses"][:]
        if plane["unscaledlens"]:
            allLenses.append(plane["unscaledlens"])

        for bf, center in allLenses:
            weights.append(1.0)
            centers.append(center/angularScale)
            ip, fp = bf.getCLParameters(angularScale, angularScale*angularScale)
            if len(ip) > 0:
                intParams += ip.tolist()
            if len(fp) > 0:
                floatParams += fp.tolist()

            print("Lens type", type(bf), "int params", ip, "float params", fp)

    intParams.append(12345) # To avoid a zero-length buffer
    floatParams.append(12345) # To avoid a zero-length buffer

    intParams = np.array(intParams, dtype=np.int32)
    floatParams= np.array(floatParams, dtype=np.float32)
    weights = np.array(weights, dtype=np.float32)
    centers = np.array(centers, dtype=np.float32)
    intParamOffsets = np.array(intParamOffsets, dtype=np.int32)
    floatParamOffsets = np.array(floatParamOffsets, dtype=np.int32)
    weightOffsets = np.array(weightOffsets, dtype=np.int32)
    planeWeights = np.array(planeWeights, dtype=np.float32)

    return allNumPlanes, DsrcAll, Dmatrix, intParams, floatParams, weights, centers, intParamOffsets, floatParamOffsets, weightOffsets, planeWeights

class CPUMultiPlaneTracer(object):
    def __init__(self, zss, zds, lensplanes, cosm):
        self.planeParams = []
        self.allWeights = []
        self.planeWeights = []
        self.hasUnscaled = []
        self.zss = zss
        self.zds = zds
        self.cosm = cosm

        for plane in lensplanes:
            params = [ {"lens": bf, "x": center[0], "y": center[1], "factor": 1, "angle": 0 } for bf, center in plane["scaledlenses"] ]
            self.allWeights.append(1.0)

            if plane["unscaledlens"]:
                bf, center = plane["unscaledlens"]
                params.append({"lens": bf, "x": center[0], "y": center[1], "factor": 1, "angle": 0 })
                self.hasUnscaled.append(True)
            else:
                self.hasUnscaled.append(False)

            for p in params:
                self.allWeights.append(1.0)

            self.planeParams.append(params)
            self.planeWeights.append(1.0)

    def getAllWeights(self):
        return self.allWeights

    def getPlaneWeights(self):
        return self.planeWeights

    def setAllWeights(self, w):
        self.allWeights = w

    def setPlaneWeights(self, w):
        self.planeWeights = w

    def trace(self, thetasArcsec):

        weightIdx = [ 0 ]
        def getNextWeight():
            w = self.allWeights[weightIdx[0]]
            weightIdx[0] += 1
            return w

        compLenses = []

        for planeIdx in range(len(self.planeParams)):
            params = self.planeParams[planeIdx]
            endIdx = len(params)-1 if self.hasUnscaled[planeIdx] else len(params)

            for paramIdx in range(endIdx):
                params[paramIdx]["factor"] = getNextWeight()*self.planeWeights[planeIdx]

            if self.hasUnscaled[planeIdx]:
                params[endIdx]["factor"] = getNextWeight()

            compLenses.append(lenses.CompositeLens(params[0]["lens"].getLensDistance(), params))

        lp = multiplane.MultiLensPlane(list(zip(compLenses, self.zds)), V(-30,-30)*ANGLE_ARCSEC, V(30,30)*ANGLE_ARCSEC, 4, 4, cosmology=self.cosm)
        betas = []
        for theta, zs in zip(thetasArcsec, self.zss):

            ip = multiplane.MultiImagePlane(lp, zs)
            beta = ip.traceTheta(theta.astype(np.double)*ANGLE_ARCSEC)/ANGLE_ARCSEC
            betas.append(beta)

        return np.array(betas)

def main2():

    cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
    zd1, zd2, zd3 = 0.3, 0.7, 1.1

    plane1 = {
        "scaledlenses": [ (lenses.PlummerLens(cosm.getAngularDiameterDistance(zd1), { "mass": 1e14*MASS_SUN, "width": 2.0*ANGLE_ARCSEC }), V(1,0)*ANGLE_ARCSEC) ],
        "unscaledlens": (lenses.MassSheetLens(cosm.getAngularDiameterDistance(zd1), { "density": 1.5 }), V(0,0)),
    }
    plane2 = {
        "scaledlenses": [ (lenses.PlummerLens(cosm.getAngularDiameterDistance(zd2), { "mass": 2e14*MASS_SUN, "width": 1.5*ANGLE_ARCSEC }), V(0,1)*ANGLE_ARCSEC) ],
        "unscaledlens": (lenses.MassSheetLens(cosm.getAngularDiameterDistance(zd2), { "density": 1.5 }), V(0,0)),
    }
    plane3 = {
        "scaledlenses": [ (lenses.CompositeLens(cosm.getAngularDiameterDistance(zd3), [
              { "lens": lenses.PlummerLens(cosm.getAngularDiameterDistance(zd3), { "mass": 1.5e14*MASS_SUN, "width": 3.5*ANGLE_ARCSEC }),
                "factor": 0.5, "angle": 0, "x": ANGLE_ARCSEC, "y": 0 },
              ]),
          V(-1,-1)*ANGLE_ARCSEC),
          (lenses.CompositeLens(cosm.getAngularDiameterDistance(zd3), [
              { "lens": lenses.PlummerLens(cosm.getAngularDiameterDistance(zd3), { "mass": 1.5e14*MASS_SUN, "width": 3.5*ANGLE_ARCSEC }),
                "factor": 0.5, "angle": 0, "x": -ANGLE_ARCSEC, "y": 0 },
              ]),
          V(-1,-1)*ANGLE_ARCSEC),
          ],
        "unscaledlens": (lenses.MassSheetLens(cosm.getAngularDiameterDistance(zd3), { "density": 1.5 }), V(0,0)),
    }
    zds, lensplanes = [zd1, zd2, zd3 ], [ plane1, plane2, plane3 ]
    code = getMultiPlaneOCLProgram(lensplanes)

    angularScale = ANGLE_ARCSEC
    thetas = (np.array([ V(5,5)*ANGLE_ARCSEC, V(-4,4)*ANGLE_ARCSEC, V(10,0)*ANGLE_ARCSEC ])/angularScale).astype(np.float32)
    zss = [ 0.5, 0.95, 2.0 ]

    #thetas = (np.array([ V(-4,4)*ANGLE_ARCSEC ])/angularScale).astype(np.float32)
    #zss = [ 0.95 ]

    betas = np.zeros(thetas.shape, dtype=np.float32)
    assert(len(thetas) == len(zss))

    allNumPlanes, DsrcAll, Dmatrix, intParams, floatParams, weights, centers, intParamOffsets, floatParamOffsets, weightOffsets, planeWeights = getOpenCLData(cosm, zss, zds, lensplanes, angularScale)

    pprint.pprint(thetas)
    pprint.pprint(betas)
    print("allNumPlanes")
    pprint.pprint(allNumPlanes)
    pprint.pprint(DsrcAll)
    pprint.pprint(Dmatrix)
    print("Int params:")
    pprint.pprint(intParams)
    print("Float params:")
    pprint.pprint(floatParams)

    print("Weights:")
    weights = np.random.random(weights.shape).astype(np.float32)
    pprint.pprint(weights)
    
    pprint.pprint(centers)
    pprint.pprint(intParamOffsets)
    pprint.pprint(floatParamOffsets)
    pprint.pprint(weightOffsets)
    
    planeWeights = np.random.random(planeWeights.shape).astype(np.float32)
    pprint.pprint(planeWeights)

    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)

    mf = cl.mem_flags
    d_thetas = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=thetas)
    d_betas = cl.Buffer(ctx, mf.WRITE_ONLY, betas.nbytes)
    d_numPlanes = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=allNumPlanes)
    d_DsrcAll = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=DsrcAll)
    d_Dmatrix = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=Dmatrix)
    d_intParams = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=intParams)
    d_floatParams = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=floatParams)
    d_weights = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=weights)
    d_centers = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=centers)
    d_intParamOffsets = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=intParamOffsets)
    d_floatParamOffsets = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=floatParamOffsets)
    d_weightOffsets = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=weightOffsets)
    d_planeWeights = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=planeWeights)

    #print(code)
    
    #print("Using code from tst.cl")
    #code = open("tst.cl", "rt").read()

    prg = cl.Program(ctx, code).build()
    calcBetas = prg.calculateBetas

    calcBetas(queue, (len(zss),), None,
             np.int32(len(zss)), d_thetas, d_betas, d_numPlanes, 
             d_DsrcAll, d_Dmatrix, d_intParams, d_floatParams, d_weights, d_centers,
             d_intParamOffsets, d_floatParamOffsets, d_weightOffsets, d_planeWeights)

    cl.enqueue_copy(queue, betas, d_betas)

    print("GPU betas:")
    pprint.pprint(betas)

    tracer = CPUMultiPlaneTracer(zss, zds, lensplanes, cosm)
    tracer.setAllWeights(weights.tolist())
    tracer.setPlaneWeights(planeWeights.tolist())
    betas = tracer.trace(thetas)
    print("CPU:")
    print(betas)

if __name__ == "__main__":
    main2()
