from grale.all import *
import pprint

def createThetaGrid(bottomLeft, topRight, numX, numY):
    thetas = np.empty([numY,numX,2], dtype=np.double)
    thetas[:,:,0], thetas[:,:,1] = np.meshgrid(np.linspace(bottomLeft[0], topRight[0], numX), 
                                               np.linspace(bottomLeft[1], topRight[1], numY))
    return thetas

def check(x):
    if not x:
        raise Exception("Not equal!")
    return True

l = lenses.PlummerLens(1000*DIST_MPC, {"mass": 1e14*MASS_SUN, "width": 2*ANGLE_ARCSEC})
bl = -V(10,10)*ANGLE_ARCSEC
tr = -bl

Nx = 511
Ny = 511

Dds = 1
Ds = 1

for threads in [ 0, 1, 2, 3, 4, 5, 6, 7, 8]:
    lenses.setNumberOfCalculationThreads(threads)
    print("Number of calculation threads:", lenses.getNumberOfCalculationThreads())
    for Nx in [ 1, 3, 17, 53, 48, 127 ]:
        for Ny in [ 1, 3, 17, 53, 48, 127 ]:

            thetas = createThetaGrid(bl, tr, Nx, Ny)

            betas = []
            for i in [ False, True ]:
                lenses.experimentalThreads = i
                betas.append(l.traceTheta(Ds, Dds, thetas))

            print("traceTheta", check(np.array_equal(betas[0], betas[1])))
            print()

            alphas = []
            for i in [ False, True ]:
                lenses.experimentalThreads = i
                alphas.append(l.getAlphaVector(thetas))

            print("getAlphaVector", check(np.array_equal(alphas[0], alphas[1])))
            print()

            derivs = []
            for i in [ False, True ]:
                lenses.experimentalThreads = i
                derivs.append(l.getAlphaVectorDerivatives(thetas))

            print("getAlphaVectorDerivatives", check(np.array_equal(derivs[0], derivs[1])))
            print()

            potentials = []
            for i in [ False, True ]:
                lenses.experimentalThreads = i
                potentials.append(l.getProjectedPotential(Ds, Dds, thetas))

            print("getProjectedPotential", check(np.array_equal(potentials[0], potentials[1])))
            print()

            denses = []
            for i in [ False, True ]:
                lenses.experimentalThreads = i
                denses.append(l.getSurfaceMassDensity(thetas))

            print("getSurfaceMassDensity", check(np.array_equal(denses[0], denses[1])))
            print()

            invmags = []
            for i in [ False, True ]:
                lenses.experimentalThreads = i
                invmags.append(l.getInverseMagnification(Ds, Dds, thetas))

            print("getInverseMagnification", check(np.array_equal(invmags[0], invmags[1])))
            print()


