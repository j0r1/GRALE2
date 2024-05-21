from grale.all import *
import random
import pickle

feedback.setDefaultFeedback("none")

# Write the random number generator state, for reproducibility
print("RNG State:")
print(random.getstate())

# Lens redshifts and cosmological model
z_lens = 0.099
z_s = 1.241

cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
cosmology.setDefaultCosmology(cosm)

# Helper function to read the input file. The images are those from
# Massey et al (https://ui.adsabs.harvard.edu/abs/2018MNRAS.477..669M/abstract)

def readInputFile():
    # This is what I used in the overlay as well
    ctrRa, ctrDec = 330.472325, -59.9454361111111

    def la(line):
        # e.g.: Ah4 330.47372 -59.94599
        idStr, ra, dec = line.split()
        ra, dec = float(ra), float(dec)
        imgNr = int(idStr[2:])
        group = idStr[1:2]
        srcNr = "A" + group # for point images, we're going to use each image group as a separate source

        x, y = images.centerOnPosition(V(ra,dec)*ANGLE_DEGREE, V(ctrRa, ctrDec)*ANGLE_DEGREE)
        return { "x": x, "y": y, "z": z_s, "srcnr": srcNr, "imgnr": imgNr }

    return images.readInputImagesFile("a3827_images.txt", True, la)

imgList = readInputFile()

regionSize = 30*ANGLE_ARCSEC
regionCenter = V(-2,0)*ANGLE_ARCSEC

# For this inversion, we'll use the input images as well as a null space

nullSize = 60*ANGLE_ARCSEC
nullSize2 = V(nullSize, nullSize)/2
null = images.createGridTriangles(regionCenter-nullSize2, regionCenter+nullSize2, 64, 64)

# Create an inversion workspace and add the data
iws = inversion.InversionWorkSpace(z_lens, regionSize=regionSize, regionCenter=regionCenter)
for i in imgList:
    iws.addImageDataToList(i["imgdata"], z_s, "pointimages")
    iws.addImageDataToList(null, z_s, "pointnullgrid")

# We're going to add two kinds of basis functions. The first is based on a regular grid
# covering the strong lensing area, and by default Plummer basis functions will be used
iws.setUniformGrid(15)
iws.addBasisFunctionsBasedOnCurrentGrid()

# Then we'll load the list of basis functions that was created in the notebook, basically
# a list of SIE models, one for each galaxy we identified. Each basis function is already
# offset to the correct center location, so in the function call below the center should
# not be set again.
sieBasisFunctions = pickle.load(open("siebasisfunctions.pickle", "rb"))
iws.addBasisFunctions([ { "lens": b, "center": V(0,0) } for b in sieBasisFunctions ])

# For the inversion we want to use all of these basis functions, so the function to
# use is 'invertBasisFunctions' (and not 'invert', which would just use basis functions
# based on the grid)
lens, fitness, fitdesc = iws.invertBasisFunctions(128)
lens.save(f"inv_hybrid.lensdata")

print("Received result with fitness", fitness, "for order:", fitdesc)

