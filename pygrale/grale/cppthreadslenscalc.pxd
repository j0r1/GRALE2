from libcpp cimport bool as cbool
from libcpp.vector cimport vector
from libcpp.string cimport string
cimport grale.vector2d as vector2d
cimport grale.gravitationallens as gravitationallens

cdef extern from "threadslenscalc.h":

    cbool threadsTraceTheta(gravitationallens.GravitationalLens &l, string &errStr,
        double Ds, double Dds,
        const double *thetaX, const double *thetaY, size_t thetaStride,
        double *betaX, double *betaY, size_t betaStride,
        size_t numElements,
        size_t nThreads)

    cbool threadsGetAlphaVector(gravitationallens.GravitationalLens &l, string &errStr,
                           const double *thetaX, const double *thetaY, size_t thetaStride,
                           double *alphaX, double *alphaY, size_t alphaStride,
                           size_t numElements,
                           size_t nThreads)
            
    cbool threadsGetAlphaVectorDerivatives(gravitationallens.GravitationalLens &l, string &errStr,
                           const double *thetaX, const double *thetaY, size_t thetaStride,
                           double *alphaXX, double *alphaYY, double *alphaXY, size_t alphaStride,
                           size_t numElements,
                           size_t nThreads)

    cbool threadsGetAlphaVectorSecondDerivatives(gravitationallens.GravitationalLens &l, string &errStr,
                           const double *thetaX, const double *thetaY, size_t thetaStride,
                           double *alphaXXX, double *alphaYYY, double *alphaXXY, double *alphaYYX,
                           size_t alphaStride,
                           size_t numElements,
                           size_t nThreads)

    cbool threadsGetInverseMagnification(gravitationallens.GravitationalLens &l, string &errStr,
                           double Ds, double Dds,
                           const double *thetaX, const double *thetaY, size_t thetaStride,
                           double *invMag, size_t invMagStride,
                           size_t numElements,
                           size_t nThreads)

    cbool threadsGetSurfaceMassDensity(gravitationallens.GravitationalLens &l, string &errStr,
                           const double *thetaX, const double *thetaY, size_t thetaStride,
                           double *dens, size_t densStride,
                           size_t numElements,
                           size_t nThreads)

    cbool threadsGetProjectedPotential(gravitationallens.GravitationalLens &l, string &errStr,
                           double Ds, double Dds,
                           const double *thetaX, const double *thetaY, size_t thetaStride,
                           double *potential, size_t potentialStride,
                           size_t numElements,
                           size_t nThreads)


