from libcpp cimport bool
from libcpp.vector cimport vector

cimport grale.vector2d
cimport grale.errut
cimport grale.polygon2d

ctypedef grale.vector2d.Vector2Dd Vector2Dd
ctypedef grale.errut.ErrorBase ErrorBase
ctypedef grale.polygon2d.Polygon2Dd Polygon2Dd

cdef extern from "grale/sourceimage.h" namespace "grale::SourceImage":
    cdef enum SourceType:
        Circle,
        Ellipse,
        Polygon,
        Discrete,
        Point

cdef extern from "grale/sourceimage.h" namespace "grale":

    cdef cppclass SourceImage(ErrorBase):

        SourceType getSourceType()
        SourceImage *createCopy()

        double getIntensity(Vector2Dd beta)
        bool isSourceInRange(Vector2Dd beta, double radius)
        double getMaxRadius()

        Vector2Dd getAngularPosition()
        void addToAngularPosition(Vector2Dd p)
        void setAngularPosition(Vector2Dd p)
        void addToAngle(double ang)	
        void setAngle(double ang)
        double getAngle()

cdef extern from "grale/circularsource.h" namespace "grale":
    cdef cppclass CircularSource(SourceImage):
        CircularSource(Vector2Dd pos, double radius, double brightnessScale)
        double getAngularRadius()
        void setAngularRadius(double r)
        void setFade(bool f)
        bool getFade()

cdef extern from "grale/ellipticalsource.h" namespace "grale":
    cdef cppclass EllipticalSource(SourceImage):
        EllipticalSource(Vector2Dd pos, double angularAxis, double eccentricty, double rotangle, double brightnessScale)
        void setFade(bool f)
        bool getFade()

cdef extern from "grale/polygonsource.h" namespace "grale":
    cdef cppclass PolygonSource(SourceImage):
        PolygonSource(Vector2Dd pos, const Polygon2Dd &polygon, double rotAngle, double brightnessScale)

cdef extern from "grale/discretesource.h" namespace "grale":
    cdef cppclass DiscreteSource(SourceImage):
        DiscreteSource(Vector2Dd angularpos, double angle, double brightnessScale)
        bool setData(const vector[double] &values, int numX, int numY, double angularWidth, double angularHeight)

cdef extern from "grale/pointsource.h" namespace "grale":
    cdef cppclass PointSource(SourceImage):
        PointSource(Vector2Dd angularpos, double brightnessScale)

