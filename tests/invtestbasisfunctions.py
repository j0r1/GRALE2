# For reproducibilty
import os
import random
os.environ["GRALE_DEBUG_SEED"] = "12345"
random.seed(12345)

from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.plotutil as plotutil
import grale.feedback as feedback
import grale.lenses as lenses
from grale.cosmology import Cosmology
import grale.images as images
import numpy as np
import matplotlib.pyplot as plt
import sys

popSize = 32
maxGen = 50
showPlots = False

feedback.setDefaultFeedback("none")
#inversion.setDefaultInverter("mpics")

pointData = """
-24.72085695	-87.40427296	1.8379
27.49850916	14.17810948
24.20549548	29.67064975

28.42028525	-52.23780402	2.1302
21.1886235	-23.67879895
1.196384135	-13.94864117
-16.78717681	-3.972760312
41.53065706	71.9816807

50.20372644	-47.07064435	2.8310
-53.33057409	-42.07215863
22.13035426	-19.29616963
-61.56069937	18.58777874
28.51342392	71.02551548

-7.309501299	-76.20927696	1.36869
23.88425492	12.48020577
12.98410719	15.79977737

-38.07680372	-32.27629926	3.1051
-18.56056319	-4.632942462
59.87129173	75.69073216

76.63073456	-39.22479261	2.1376
-20.93129261	3.389082883
-35.27228316	17.34904039

-32.95280779	-64.40176566	2.5102
64.5943462	-35.72726102
16.54139764	-1.488663532
0.9954536365	-1.167836662
-15.77443336	0.4511665107
-46.81516718	32.04952996
39.88120213	59.4284112

-30.97837902	-73.49947513	1.7204
28.4747397	10.30246005
47.91205058	32.0326518

-64.73715927	-60.8487126	3.4344
61.88479429	5.053116483
30.9402014	7.874761056

-20.26552689	0.2918478286	2.0480
-44.9058689	10.37961941
72.08870539	44.05628017

-54.10336245	-67.95715467	2.5561
30.44231254	9.113465759
59.49758414	14.07971587

-61.28338007	-40.56482628	3.3279
56.01285194	-39.39177112
23.54727476	-17.99304283
-57.8598424	35.14643033
13.60282371	72.3803249

-20.34211625	1.175010103	2.3608
-43.01975403	14.68476623
81.45607978	30.77891962

25.19265142	-90.01278091	2.8322
-22.02145039	10.18050396
-25.42494506	16.10820935

-9.731165787	-34.78165139	2.8319
-16.39188153	-12.63667786
54.72215837	84.1282582

-21.70735217	-84.12504853	2.4716
-12.06767687	8.766316014
25.12333893	10.48491248
-3.292798587	10.93851762
-19.67473486	35.76443409

-71.96458367	-61.98898097	3.1210
34.20162797	8.778215329
50.7348093	9.249756314

-21.03927634	-0.8220949323	1.9547
-43.46606829	1.318453962
72.03580578	51.03686822

8.921826202	-71.16791105	2.9976
30.0613313	-68.38836807
-18.60631484	-0.2792288361
-52.5929001	14.90107726
64.73092626	55.63898176

25.52338997	-85.74702361	2.3711
-21.81376938	8.271518002
-26.94970019	16.35484607

53.56965743	-43.54761517	3.4573
-61.33258478	-37.6570158
23.15899998	-19.25245575
-62.46838245	26.89455438
19.53228144	73.75671277

47.28119499	-51.72166526	2.7501
-43.55867768	-49.49272263
19.30124747	-16.9427861
-0.4142903709	-6.376420558
-15.50428292	-1.884038019
-60.06804983	14.08154015
37.96515147	68.83852853

-64.12416086	-82.79545034	3.2336
36.2571622	13.72493334
39.32244461	15.59162438
"""

z_lens = 0.4

#for useWeights in [ True ]:
for useWeights in [ False ]:
    iws = inversion.InversionWorkSpace(z_lens, 250*ANGLE_ARCSEC, cosmology=Cosmology(0.7, 0.27, 0, 0.73))

    imgList = images.readInputImagesFile(pointData, True) 
    for i in imgList:
        iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")

    iws.setUniformGrid(15)
    if showPlots:
        plotutil.plotSubdivisionGrid(iws.getGrid())
        plt.show()
    lens1, fitness1, fitdesc = iws.invert(popSize, maximumGenerations=maxGen, rescaleBasisFunctions=useWeights)

    iws.setSubdivisionGrid(lens1, 300, 400)
    if showPlots:
        plotutil.plotSubdivisionGrid(iws.getGrid())
        plt.show()

    gridMassEst = iws.estimateStrongLensingMass()

    # This lens model function creates the same basis functions as the grid based
    # approach uses
    def lensModelFunction(opType, opInfo, parameters = None):
        if opType == "start":
            numCells = len(opInfo["grid"])
            cellMass = gridMassEst/numCells
            return { "cellmass": cellMass }

        assert(opType == "add")

        factor = 1.7
        width = opInfo["size"]*factor
        mass = parameters["cellmass"]
        #mass *= (width/ANGLE_ARCSEC)**2 # smaller cells have less mass

        lens = lenses.PlummerLens(iws.getLensDistance(), { "width": width, "mass": mass })
        return lens, mass


    iws.clearBasisFunctions()
    iws.addBasisFunctionsBasedOnCurrentGrid(lensModelFunction)

    lens2a, fitness2a, fitdesc = iws.invert(popSize, maximumGenerations=maxGen, rescaleBasisFunctions=useWeights)
    lens2, fitness2, fitdesc = iws.invertBasisFunctions(popSize, maximumGenerations=maxGen)

    for idx, fitness, lens in [ (1, fitness1, lens1), (2, fitness2, lens2), ("2a", fitness2a, lens2a) ]:
        print(f"Fitness for lens {idx}: {fitness}")
        print(f"Recalculated:", iws.calculateFitness(lens))
        print()

    print("Press enter to continue")
    input()

