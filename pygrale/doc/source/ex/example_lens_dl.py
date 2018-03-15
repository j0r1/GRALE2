import grale.lenses as lenses
import urllib.request

modelUrl = "https://github.com/j0r1/LensModels/blob/master/2008MNRAS.389..415L_CL0024/mainresult_fig2_fig3.lensdata?raw=true"
modelData = urllib.request.urlopen(modelUrl).read()
cl0024Model = lenses.GravitationalLens.fromBytes(modelData)

