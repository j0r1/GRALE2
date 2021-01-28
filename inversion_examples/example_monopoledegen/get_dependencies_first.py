import grale.lenses as lenses
import urllib.request
import os

files = [ "mainresult_fig2_fig3.lensdata",
          "source1entire.imgdata",
          "source2entire.imgdata" ]

modelUrl = "https://github.com/j0r1/LensModels/blob/master/2008MNRAS.389..415L_CL0024/{}?raw=true"

for f in files:
    if not os.path.exists(f):
        url = modelUrl.format(f)
        print("Downloading", f)

        data = urllib.request.urlopen(url).read()
        open(f, "wb").write(data)
    else:
        print("File", f, "already exists, skipping download")


if os.path.exists("approxlens.lensdata"):
    print("approxlens.lensdata already exists, skipping generation")
else:
    print("Generating approxlens.lensdata")

    from grale.all import *

    l = lenses.GravitationalLens.load("mainresult_fig2_fig3.lensdata")
    li = plotutil.LensInfo(l, size=2*ANGLE_ARCMIN)
    lp = li.getLensPlane()
    approx = lp.createDeflectionGridLens()
    approx.save("approxlens.lensdata")
