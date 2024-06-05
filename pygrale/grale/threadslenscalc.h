#pragma once

#include <grale/gravitationallens.h>

bool threadsTraceTheta(grale::GravitationalLens &l, std::string &errStr,
		               double Ds, double Dds,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *betaX, double *betaY, size_t betaStride,
					   size_t numElements,
					   size_t nThreads);

bool threadsGetAlphaVector(grale::GravitationalLens &l, std::string &errStr,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *alphaX, double *alphaY, size_t alphaStride,
					   size_t numElements,
					   size_t nThreads);
		
bool threadsGetAlphaVectorDerivatives(grale::GravitationalLens &l, std::string &errStr,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *alphaXX, double *alphaYY, double *alphaXY, size_t alphaStride,
					   size_t numElements,
					   size_t nThreads);

bool threadsGetAlphaVectorSecondDerivatives(grale::GravitationalLens &l, std::string &errStr,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *alphaXXX, double *alphaYYY, double *alphaXXY, double *alphaYYX,
					   size_t alphaStride,
					   size_t numElements,
					   size_t nThreads);

bool threadsGetInverseMagnification(grale::GravitationalLens &l, std::string &errStr,
		               double Ds, double Dds,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *invMag, size_t invMagStride,
					   size_t numElements,
					   size_t nThreads);

bool threadsGetSurfaceMassDensity(grale::GravitationalLens &l, std::string &errStr,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *dens, size_t densStride,
					   size_t numElements,
					   size_t nThreads);

bool threadsGetProjectedPotential(grale::GravitationalLens &l, std::string &errStr,
		               double Ds, double Dds,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *potential, size_t potentialStride,
					   size_t numElements,
					   size_t nThreads);

