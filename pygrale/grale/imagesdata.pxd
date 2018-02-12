from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

cimport grale.vector2d
cimport grale.triangleindices

ctypedef grale.vector2d.Vector2Dd Vector2Dd
ctypedef grale.triangleindices.TriangleIndices TriangleIndices

cdef extern from "grale/imagesdata.h" namespace "grale":
    cdef cppclass ImagesData:
        ImagesData()
        string getErrorString()
        
        bool create(int numImages, bool intensities, bool shearInfo)
        int addImage();
        int addPoint(int imageNumber, Vector2Dd point)
        int addPoint(int imageNumber, Vector2Dd point, double intensity)
        int addPoint(int imageNumber, Vector2Dd point, double shearComponent1, double shearComponent2)
        int addPoint(int imageNumber, Vector2Dd point, double intensity, double shearComponent1, double shearComponent2)

        int addGroup()
        bool addGroupPoint(int groupNumber, int imageIndex, int pointIndex)
        bool addTimeDelayInfo(int imageIndex, int pointIndex, double timeDelay)
        bool addTriangle(int imageNumber, int index1, int index2, int index3)

        bool load(string fname)
        bool save(const string &fname)
        
        bool hasIntensities()
        bool hasShearInfo()
        bool hasTimeDelays()
    
        int getNumberOfImages()
        int getNumberOfImagePoints(int i)
        Vector2Dd getImagePointPosition(int image, int point)
        double getImagePointIntensity(int image, int point)
        double getShearComponent1(int image, int point)
        double getShearComponent2(int image, int point)
    
        int getNumberOfGroups()
        int getNumberOfGroupPoints(int group)
        void getGroupPointIndices(int group, int pointnr, int *img, int *point)
        
        int getNumberOfTimeDelays()
        void getTimeDelay(int index, int *pImg, int *pPoint, double *pDelay)

        bool hasTriangulation()    
        bool getTriangles(int image, vector[TriangleIndices] &triangles)
        void clearTriangulation()

        Vector2Dd getTopRightCorner()
        Vector2Dd getBottomLeftCorner()

        void centerOnPosition(double ra, double dec)
        void uncenterOnPosition(double ra, double dec)
        void subtractIntensity(double v)

