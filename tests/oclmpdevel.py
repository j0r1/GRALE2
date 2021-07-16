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

def traceToSourcePlane_0(theta, numPlanes, Dsrc, Dmatrix, getAlpha):

    T = [ None for _ in range(numPlanes) ]

    for j in range(0,numPlanes):
        T[j] = theta.copy()
        for i in range(0, j):
            alpha = getAlpha(i, T[i])
            T[j] -= Dmatrix[i+1,j]/Dmatrix[0,j] * alpha

    beta = theta.copy()
    for i in range(0, numPlanes):
        alpha = getAlpha(i, T[i])
        beta -= Dsrc[i+1]/Dsrc[0] * alpha

    return beta

def traceToSourcePlane(theta, numPlanes, Dsrc, Dmatrix, getAlpha):

    T = [ None for _ in range(numPlanes) ]
    alphas = [ None for _ in range(numPlanes-1) ]

    for j in range(0,numPlanes):
        T[j] = theta.copy()
        for i in range(0, j-1):
            T[j] -= Dmatrix[i+1,j]/Dmatrix[0,j] * alphas[i]

        i = j-1
        if i >= 0:
            alpha = getAlpha(i, T[i])
            dfrac = Dmatrix[i+1,j]/Dmatrix[0,j]
            #print(f"pyalpha[{i}] = {alpha} T[{i}] = {T[i]} dfrac = {dfrac}")
            alphas[i] = alpha
            T[j] -= dfrac * alpha

    beta = theta.copy()
    for i in range(0, numPlanes-1):
        beta -= Dsrc[i+1]/Dsrc[0] * alphas[i]

    i = numPlanes-1
    if i >= 0:
        alpha = getAlpha(i, T[i])
        dfrac = Dsrc[i+1]/Dsrc[0]
        #print(f"pyalpha[{i}] = {alpha} T[{i}] = {T[i]} dfrac = {dfrac}")
        beta -= dfrac * alpha

    return beta

def getNeededPlanes(zs, zds):
    needed = 0
    for zd in zds:
        if zs > zd:
            needed += 1
        else:
            break
    return needed

def testTrace(zs, zds, lensplanes, cosm):

    Dmatrix = getDistanceMatrix(zds, cosm)
    Dsrc, usedPlanes = getDistanceVector(zs, zds, cosm)
    pprint.pprint(Dmatrix)
    pprint.pprint(Dsrc)
    
    getAlpha = lambda idx, theta : lensplanes[idx].getAlphaVector(theta*ANGLE_ARCSEC)/ANGLE_ARCSEC

    numPlanes = getNeededPlanes(zs, zds)
    print(zs, zds, numPlanes)
    theta = V(5,7)
    beta = traceToSourcePlane(theta, numPlanes, Dsrc, Dmatrix, getAlpha)

    if numPlanes > 0:
        lp = multiplane.MultiLensPlane(list(zip(lensplanes, zds)), V(-30,-30)*ANGLE_ARCSEC, V(30,30)*ANGLE_ARCSEC, 4, 4, cosmology=cosm)
        ip = multiplane.MultiImagePlane(lp, zs)
        beta2 = ip.traceTheta(theta*ANGLE_ARCSEC)/ANGLE_ARCSEC
        print("theta = ", theta, "beta = ", beta, "beta2 = ", beta2)
    else:
        print("theta = ", theta, "beta = ", beta)


def main():
    cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
    zd1, zd2, zd3 = 0.3, 0.7, 1.1
    zds, lensplanes = [zd1, zd2, zd3 ], [ 
            lenses.PlummerLens(cosm.getAngularDiameterDistance(zd1), { "mass": 1e14*MASS_SUN, "width": 2.0*ANGLE_ARCSEC }),
            lenses.PlummerLens(cosm.getAngularDiameterDistance(zd2), { "mass": 2e14*MASS_SUN, "width": 1.5*ANGLE_ARCSEC }),
            lenses.PlummerLens(cosm.getAngularDiameterDistance(zd3), { "mass": 1.5e14*MASS_SUN, "width": 3.5*ANGLE_ARCSEC }),
    ]

    testTrace(0.1, zds, lensplanes, cosm)
    testTrace(0.5, zds, lensplanes, cosm)
    testTrace(1.0, zds, lensplanes, cosm)
    testTrace(2.5, zds, lensplanes, cosm)
    
def oclTest():

	a_np = np.random.rand(50000).astype(np.float32)
	b_np = np.random.rand(50000).astype(np.float32)

	ctx = cl.create_some_context()
	queue = cl.CommandQueue(ctx)

	mf = cl.mem_flags
	a_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=a_np)
	b_g = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_np)

	prg = cl.Program(ctx, """
	__kernel void sum(
		__global const float *a_g, __global const float *b_g, __global float *res_g)
	{
	  int gid = get_global_id(0);
	  res_g[gid] = a_g[gid] + b_g[gid];
	}
	""").build()

	res_g = cl.Buffer(ctx, mf.WRITE_ONLY, a_np.nbytes)
	knl = prg.sum  # Use this Kernel object for repeated calls
	knl(queue, a_np.shape, None, a_g, b_g, res_g)

	res_np = np.empty_like(a_np)
	cl.enqueue_copy(queue, res_np, res_g)

	# Check on CPU with Numpy:
	print(res_np - (a_np + b_np))
	print(np.linalg.norm(res_np - (a_np + b_np)))
	assert np.allclose(res_np, a_np + b_np)

# TODO: same types should be consecutive, will be much more performant
def getAlphaCodeForPlane(functionName, plane, neededBasisFunctionCode):

    code = """
float2 {functionName}(float2 theta, __global const int *pIntParams, __global const float *pFloatParams, __global const float *pWeights,
                      __global const float *pCenters)
{{
    float2 alpha = (float2)(0, 0);
    //printf("{functionName}\\n");
""".format(functionName=functionName)

    def addCodeForType(t, num):
        entry = neededBasisFunctionCode[t]
        code = """
    for (int i = 0 ; i < {Nlenses} ; i++)
    {{
        const float2 center = (float2)(pCenters[0], pCenters[1]);
        //printf("center = %g %g\\n", center.x, center.y);
        LensQuantities l = {lensFunctionName}(theta-center, pIntParams, pFloatParams);
        float w = *pWeights;
        alpha.x += w*l.alphaX;
        alpha.y += w*l.alphaY;
    
        pWeights++;
        pCenters += 2;
""".format(Nlenses=num, lensFunctionName=entry["functionname"])
        if entry["intcount"]:
            code += "        pIntParams += {};\n".format(entry["intcount"])
        if entry["floatcount"]:
            code += "        pFloatParams += {};\n".format(entry["floatcount"])
        code += """
    }
"""
        return code

    prevBfType, prevBfCounts = None, 0
    for bf,center in plane:
        t = type(bf)
        if not t in neededBasisFunctionCode:
            fnName, fnCode = bf.getCLProgram(False, False)
            iCnt, fCnt = bf.getCLParameterCounts()
            neededBasisFunctionCode[t] = {
                "functionname": fnName,
                "functioncode": fnCode,
                "intcount": iCnt,
                "floatcount": fCnt,
            }
        else:
            fnName = neededBasisFunctionCode[t]["functionname"]

        if prevBfType != t:
            if prevBfCounts:
                code += addCodeForType(prevBfType, prevBfCounts)
            prevBfType, prevBfCounts = t, 1
        else:
            prevBfCounts += 1

    if prevBfCounts:
        code += addCodeForType(prevBfType, prevBfCounts)

    code += """
    return alpha;
}""";
    return code

# TODO: separate scale parameter ? mass sheet basis function that doesn't use this scale?

def getMultiPlaneTraceCode(lensPlanes):
    maxPlanes = len(lensPlanes)

    alphaCode = ""
    neededBasisFunctionCode = { }

    code = """

#define MAXPLANES {maxPlanes}

float2 getAlpha(int lpIdx, float2 theta, __global const int *pAllIntParams, __global const float *pAllFloatParams,
                __global const float *pAllWeights, __global const float *pAllCenters, __global const int *pPlaneIntParamOffsets,
                __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets)
{{

""".format(maxPlanes=maxPlanes)
    
    for i in range(maxPlanes):
        alphaCode += getAlphaCodeForPlane("getAlpha_{}".format(i), lensPlanes[i], neededBasisFunctionCode)
        code += """
    if (lpIdx == {i})
        return getAlpha_{i}(theta, pAllIntParams + pPlaneIntParamOffsets[{i}], pAllFloatParams + pPlaneFloatParamOffsets[{i}],
                            pAllWeights + pPlaneWeightOffsets[{i}], pAllCenters + pPlaneWeightOffsets[{i}]*2);
""".format(i = i);

    code += """

    return (float2)(0.0f/0.0f, 0.0f/0.0f);
}

float2 multiPlaneTrace(float2 theta, int numPlanes, __global const float *Dsrc, __global const float *Dmatrix,
                __global const int *pAllIntParams, __global const float *pAllFloatParams, __global const float *pAllWeights,
                __global const float *pAllCenters, __global const int *pPlaneIntParamOffsets,
                __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets)
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
            const float2 alpha = getAlpha(i, T[i], pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters, pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets);
            alphas[i] = alpha;
	        const float dfrac = Dmatrix[(i+1)*MAXPLANES + j]/Dmatrix[0 + j];
            //printf("alpha[%d] = %g %g T[%d] = %g %g dfrac = %g\\n",i,alpha.x,alpha.y,i,T[i].x,T[i].y,dfrac);

            T[j] -= dfrac * alpha;
        }
    }

    float2 beta = theta;
    for (int i = 0 ; i < numPlanes-1 ; i++)
        beta -= (Dsrc[i+1]/Dsrc[0]) * alphas[i];

    const int i = numPlanes-1;
    if (i >= 0)
    {
        const float2 alpha = getAlpha(i, T[i], pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters, pPlaneIntParamOffsets, pPlaneFloatParamOffsets, pPlaneWeightOffsets);
        const float dfrac = Dsrc[i+1]/Dsrc[0];
        //printf("alpha[%d] = %g %g T[%d] = %g %g dfrac = %g\\n",i,alpha.x,alpha.y,i,T[i].x,T[i].y,dfrac);
        beta -= dfrac * alpha;
    }
    return beta;
}
"""

    allCode = ""
    for t in neededBasisFunctionCode:
        allCode += neededBasisFunctionCode[t]["functioncode"]
    allCode += alphaCode
    allCode += code

    return allCode

def getMultiPlaneOCLProgram(lensplanes):

    code = lensplanes[0][0][0].getCLLensQuantitiesStructure(False, False) + getMultiPlaneTraceCode(lensplanes) + """
__kernel void calculateBetas(const int numPoints, __global const float *pThetas, __global float *pBetas, 
                                __global const int *pNumPlanes, __global const float *DsrcAll, __global const float *Dmatrix,
                                __global const float *pAllIntParams, __global const float *pAllFloatParams,
                                __global const float *pAllWeights, __global const float *pAllCenters,
                                __global const int *pPlaneIntParamOffsets,
                                __global const int *pPlaneFloatParamOffsets, __global const int *pPlaneWeightOffsets)
{
    const int i = get_global_id(0);
    if (i >= numPoints)
        return;

    //printf("pAllCenters = %g %g / %g %g / %g %g\\n", pAllCenters[0], pAllCenters[1], pAllCenters[2], pAllCenters[3], pAllCenters[4], pAllCenters[5]);
    //printf("pPlaneIntParamOffsets = %d %d %d\\n", pPlaneIntParamOffsets[0], pPlaneIntParamOffsets[1], pPlaneIntParamOffsets[2]);
    //printf("pPlaneFloatParamOffsets = %d %d %d\\n", pPlaneFloatParamOffsets[0], pPlaneFloatParamOffsets[1], pPlaneFloatParamOffsets[2]);
    //printf("pPlaneWeightOffsets = %d %d %d\\n", pPlaneWeightOffsets[0], pPlaneWeightOffsets[1], pPlaneWeightOffsets[2]);
    //printf("pNumPlanes = %d %d %d\\n", pNumPlanes[0], pNumPlanes[1], pNumPlanes[2]);


    const float2 theta = (float2)(pThetas[i*2+0], pThetas[i*2+1]);
    // Each Dsrc is vector of MAXPLANES+1 length
    //printf("theta[%d] = %g %g numPlanes[%d] = %d\\n", i, theta.x, theta.y, i, pNumPlanes[i]);
    __global const float *Dsrc = DsrcAll + (MAXPLANES+1)*i;
    //printf("Dsrc[%d] = %g %g %g %g\\n", i, Dsrc[0], Dsrc[1], Dsrc[2], Dsrc[3]);
    const float2 beta = multiPlaneTrace(theta, pNumPlanes[i], Dsrc, Dmatrix,
                                        pAllIntParams, pAllFloatParams, pAllWeights, pAllCenters,
                                        pPlaneIntParamOffsets,
                                        pPlaneFloatParamOffsets, pPlaneWeightOffsets);

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

    for plane in lensplanes:
        intParamOffsets.append(len(intParams))
        floatParamOffsets.append(len(floatParams))
        weightOffsets.append(len(weights))

        for bf, center in plane:
            weights.append(1.0)
            centers.append(center/angularScale)
            ip, fp = bf.getCLParameters(angularScale, angularScale*angularScale)
            if len(ip) > 0:
                intParams += ip.tolist()
            if len(fp) > 0:
                floatParams += fp.tolist()

    intParams.append(0) # To avoid a zero-length buffer
    floatParams.append(0) # To avoid a zero-length buffer

    intParams = np.array(intParams, dtype=np.float32)
    floatParams= np.array(floatParams, dtype=np.float32)
    weights = np.array(weights, dtype=np.float32)
    centers = np.array(centers, dtype=np.float32)
    intParamOffsets = np.array(intParamOffsets, dtype=np.int32)
    floatParamOffsets = np.array(floatParamOffsets, dtype=np.int32)
    weightOffsets = np.array(weightOffsets, dtype=np.int32)

    return allNumPlanes, DsrcAll, Dmatrix, intParams, floatParams, weights, centers, intParamOffsets, floatParamOffsets, weightOffsets

def traceCPU(thetas, zss, zds, lensplanes, cosm):

    print("CPU:")

    compLenses = []
    for plane in lensplanes:
        params = [ {"lens": bf, "x": center[0], "y": center[1], "factor": 1, "angle": 0 } for bf, center in plane ]
        compLenses.append(lenses.CompositeLens(params[0]["lens"].getLensDistance(), params))

    lp = multiplane.MultiLensPlane(list(zip(compLenses, zds)), V(-30,-30)*ANGLE_ARCSEC, V(30,30)*ANGLE_ARCSEC, 4, 4, cosmology=cosm)
    for theta, zs in zip(thetas, zss):

        ip = multiplane.MultiImagePlane(lp, zs)
        beta = ip.traceTheta(theta.astype(np.double)*ANGLE_ARCSEC)

        beta /= ANGLE_ARCSEC

        print("theta {} {} => beta {} {}".format(theta[0],theta[1],beta[0],beta[1]))

        Dmatrix = getDistanceMatrix(zds, cosm)
        Dsrc, usedPlanes = getDistanceVector(zs, zds, cosm)
        
        getAlpha = lambda idx, theta : compLenses[idx].getAlphaVector(theta.astype(np.double)*ANGLE_ARCSEC)/ANGLE_ARCSEC

        numPlanes = getNeededPlanes(zs, zds)
        beta = traceToSourcePlane(theta, numPlanes, Dsrc, Dmatrix, getAlpha)
        print("theta2 {} {} => beta2 {} {}".format(theta[0],theta[1],beta[0],beta[1]))

def main2():

    cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
    zd1, zd2, zd3 = 0.3, 0.7, 1.1
    zds, lensplanes = [zd1, zd2, zd3 ], [ 
        [ (lenses.PlummerLens(cosm.getAngularDiameterDistance(zd1), { "mass": 1e14*MASS_SUN, "width": 2.0*ANGLE_ARCSEC }), V(1,0)*ANGLE_ARCSEC) ],
        [ (lenses.PlummerLens(cosm.getAngularDiameterDistance(zd2), { "mass": 2e14*MASS_SUN, "width": 1.5*ANGLE_ARCSEC }), V(0,1)*ANGLE_ARCSEC) ],
        [ (lenses.PlummerLens(cosm.getAngularDiameterDistance(zd3), { "mass": 1.5e14*MASS_SUN, "width": 3.5*ANGLE_ARCSEC }), V(-1,-1)*ANGLE_ARCSEC) ],
    ]

    code = getMultiPlaneOCLProgram(lensplanes)
    #    [ 
    #        lenses.PlummerLens(1000*DIST_MPC, { "mass": 1e14*MASS_SUN, "width": 1*ANGLE_ARCSEC }),
    #        lenses.PlummerLens(1000*DIST_MPC, { "mass": 0.5e14*MASS_SUN, "width": 2*ANGLE_ARCSEC }),
    #        lenses.SISLens(1000*DIST_MPC, { "velocityDispersion": 300000 }),
    #        lenses.MassSheetLens(1000*DIST_MPC, { "density": 3 }),
    #    ],
    #    [ 
    #        lenses.PlummerLens(1200*DIST_MPC, { "mass": 1e14*MASS_SUN, "width": 1*ANGLE_ARCSEC }),
    #        lenses.SISLens(1200*DIST_MPC, { "velocityDispersion": 300000 }),
    #    ]
    #]
    #)

    angularScale = ANGLE_ARCSEC
    thetas = (np.array([ V(5,5)*ANGLE_ARCSEC, V(-4,4)*ANGLE_ARCSEC, V(10,0)*ANGLE_ARCSEC ])/angularScale).astype(np.float32)
    zss = [ 0.5, 0.95, 2.0 ]

    #thetas = (np.array([ V(-4,4)*ANGLE_ARCSEC ])/angularScale).astype(np.float32)
    #zss = [ 0.95 ]

    betas = np.zeros(thetas.shape, dtype=np.float32)
    assert(len(thetas) == len(zss))

    allNumPlanes, DsrcAll, Dmatrix, intParams, floatParams, weights, centers, intParamOffsets, floatParamOffsets, weightOffsets = getOpenCLData(cosm, zss, zds, lensplanes, angularScale)

    pprint.pprint(thetas)
    pprint.pprint(betas)
    print("allNumPlanes")
    pprint.pprint(allNumPlanes)
    pprint.pprint(DsrcAll)
    pprint.pprint(Dmatrix)
    pprint.pprint(intParams)
    pprint.pprint(floatParams)
    pprint.pprint(weights)
    pprint.pprint(centers)
    pprint.pprint(intParamOffsets)
    pprint.pprint(floatParamOffsets)
    pprint.pprint(weightOffsets)

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

    #print(code)
    prg = cl.Program(ctx, code).build()
    calcBetas = prg.calculateBetas

    calcBetas(queue, (len(zss),), None,
             np.int32(len(zss)), d_thetas, d_betas, d_numPlanes, 
             d_DsrcAll, d_Dmatrix, d_intParams, d_floatParams, d_weights, d_centers,
             d_intParamOffsets, d_floatParamOffsets, d_weightOffsets)

    cl.enqueue_copy(queue, betas, d_betas)

    print("Resulting betas:")
    pprint.pprint(betas)

    traceCPU(thetas, zss, zds, lensplanes, cosm)
    return





if __name__ == "__main__":
    main2()
