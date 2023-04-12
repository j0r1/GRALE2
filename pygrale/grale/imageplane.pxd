from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

cimport grale.vector2d
cimport grale.sourceimage
cimport grale.lensplane
cimport grale.gravitationallens
ctypedef grale.vector2d.Vector2Dd Vector2Dd
ctypedef grale.sourceimage.SourceImage SourceImage
ctypedef grale.lensplane.LensPlane LensPlane
ctypedef grale.gravitationallens.GravitationalLens GravitationalLens

cdef extern from "grale/imageplane.h" namespace "grale":
    cdef cppclass ImagePlane:
        ImagePlane()
        string getErrorString()
        bool init(const LensPlane *lensplane, double Ds, double Dds)
        bool isInit()
        double getDs()
        double getDds()

        int getNumXPoints() const
        int getNumYPoints() const
        int getNumXPixels() const
        int getNumYPixels() const

        Vector2Dd getBottomLeft() const	
        Vector2Dd getTopRight() const
        double getXStep() const
        double getYStep() const

        double getImageIntensityAccurate(vector[SourceImage *] &s, int x, int y, int subSamples) 
        double getSourceIntensityAccurate(vector[SourceImage *] &s, int x, int y, int subSamples)

        vector[vector[Vector2Dd]] &getCriticalLineSegments()
        vector[vector[Vector2Dd]] &getCausticSegments()

        bool traceBeta(Vector2Dd beta, vector[Vector2Dd] &thetaPoints)
        bool traceThetaApproximately(Vector2Dd theta, Vector2Dd *pBeta)

        @staticmethod
        bool staticTraceBetaApproximately(Vector2Dd beta, vector[Vector2Dd] &thetaPoints, vector[Vector2Dd] &betaMap, Vector2Dd bottomLeft, Vector2Dd topRight, int numX, int numY, string &errorString)
        @staticmethod
        bool staticRefinePosition(const GravitationalLens &lens, double Ds, double Dds, Vector2Dd beta, Vector2Dd startTheta, Vector2Dd &theta, int numIterations, string &errStr)
