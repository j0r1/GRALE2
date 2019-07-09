from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cimport grale.vector2d as vector2d
cimport grale.errut as errut

ctypedef vector2d.Vector2Df Vector2Df

cdef extern from "grale/multiplanecuda.h" namespace "grale::MultiPlaneCUDA":
    cdef cppclass PlummerInfo:
        PlummerInfo(Vector2Df position, float widthInAngularUnit, double initialMass)
        Vector2Df m_position
        float m_widthInAngularUnit
        double m_initialMass

cdef extern from "grale/multiplanecuda.h" namespace "grale":
   
    cdef cppclass MultiPlaneCUDA(errut.ErrorBase):
        bool init(string &libraryPath,
            double angularUnit,
            double h, double W_m, double W_r, double W_v, double w,
            vector[float] &lensRedshifts,
            vector[vector[PlummerInfo]] &fixedPlummerParameters, 
            vector[float] &sourceRedshifts,
            vector[vector[Vector2Df]] &theta)

        bool calculateSourcePositions(vector[vector[float]] &massFactors, vector[float] &sheetDensities)
        const vector[Vector2Df] *getSourcePositions(int srcIdx)

