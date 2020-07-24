import grale.lenses as lenses
import grale.multiplane as multiplane
import grale.cosmology as cosmology
from grale.constants import *
import numpy as np

V = lambda x, y: np.array([x, y], dtype=np.double)

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
D = cosm.getAngularDiameterDistance

cosmology.setDefaultCosmology(cosm)
z1, z2 = 0.5, 1.5

sourceRedshifts = [ 1.25, 2.0 ]
allThetas = [ [ V(10,8)*ANGLE_ARCSEC ] , [ V(-12,2)*ANGLE_ARCSEC ] ]

for factor in [ 0.5, 1.0, 2.0 ]:
    print("Using factor", factor)
    lensesAndRedshifts = [
        (lenses.CompositeLens(D(z1), [
            { "lens": lenses.PlummerLens(D(z1), { "mass": 1e14*MASS_SUN, "width": 5*ANGLE_ARCSEC }),
              "x": 0, "y": 0, "angle": 0, "factor": factor }
            ]), z1),
        (lenses.CompositeLens(D(z2), [
            { "lens": lenses.PlummerLens(D(z2), { "mass": 5e13*MASS_SUN, "width": 3*ANGLE_ARCSEC }),
              "x": 1*ANGLE_ARCSEC, "y": -1*ANGLE_ARCSEC, "angle": 0, "factor": factor }
            ]), z2),
        ]

    lensPlane = multiplane.MultiLensPlane(lensesAndRedshifts, V(-30,-30)*ANGLE_ARCSEC, V(30, 30)*ANGLE_ARCSEC, 511, 511)

    for idx, zs, thetas in zip(range(len(sourceRedshifts)), sourceRedshifts, allThetas):
        print("Positions for source", idx)

        imgPlane = multiplane.MultiImagePlane(lensPlane, zs)
        betas = imgPlane.traceTheta(np.array(thetas))
        for bx, by in betas/ANGLE_ARCSEC:
            print(f"{bx:.6f} {by:.6f}")
