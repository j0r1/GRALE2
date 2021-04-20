from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr

cimport grale.vector2d as vector2d
cimport grale.serut as serut
ctypedef vector2d.Vector2Dd Vector2Dd

cdef extern from "grale/lensplane.h" namespace "grale":
    cdef cppclass LensPlane:
        LensPlane()
        string getErrorString()

        bool init(serut.SerializationInterface &lensData, Vector2Dd bottomLeft, Vector2Dd topRight, int numx, int numy)
        bool init(serut.SerializationInterface &lensData, Vector2Dd bottomLeft, Vector2Dd topRight, int numx, int numy, serut.SerializationInterface &renderedData)
        bool isInit()

        bool scaleDeflections(double factor)
        int getNumXPoints() const
        int getNumYPoints() const
        Vector2Dd getBottomLeft() const    
        Vector2Dd getTopRight() const
        double getXStep() const
        double getYStep() const

        Vector2Dd getAlpha(int xpos, int ypos) const
        void getAlphaDerivatives(int xpos, int ypos, double &axx, double &ayy, double &axy) const

        Vector2Dd getIndexCoordinate(int xpos,int ypos) const

        @staticmethod
        bool load(const string &fname, unique_ptr[LensPlane] &ip, string &errstr)
        bool save(const string &fname) const
        @staticmethod
        bool read(serut.SerializationInterface &si, unique_ptr[LensPlane] &ip, string &errstr)
        bool write(serut.SerializationInterface &si) const

cdef extern from "pylensplane.h":
    cdef cppclass PyLensPlane(LensPlane):
        PyLensPlane(object o)

        @staticmethod
        bool createDeflectionGridLens(LensPlane *lp, vector[unsigned char] &b, string &errStr)
        @staticmethod
        bool getLensBytes(LensPlane *lp, vector[unsigned char] &b, string &errStr)


