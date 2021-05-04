#include "threadslenscalc.h"
#include <thread>
#include <vector>
#include <cassert>

using namespace std;
using namespace grale;

vector<size_t> getElementsPerThread(size_t nThreads, size_t numElements)
{
	vector<size_t> elementsInThread(nThreads, 0);

	size_t elemPerThread = 1;
	if (numElements % nThreads == 0) // nicely divisible
		elemPerThread = numElements/nThreads;
	else
		elemPerThread = numElements/nThreads + 1;
	
	size_t elemLeft = numElements;
	for (auto &e : elementsInThread)
	{
		if (elemPerThread >= elemLeft)
		{
			e = elemLeft;
			elemLeft = 0;
			break;
		}
		else
		{
			e = elemPerThread;
			elemLeft -= elemPerThread;
		}
	}
	assert(elemLeft == 0);
	return elementsInThread;
}

template<class X, class C, class U>
bool runThreads(GravitationalLens &l, string &errStr,
					   size_t numElements,
					   size_t nThreads,
					   X startSingleThread, C createOneThread, U updatePointers)
{
	if (numElements == 0)
		return true; // nothing to do

	auto E = [&errStr](const string &s)
	{
		errStr = s;
		return false;
	};

	if (nThreads < 1)
		return E("No threads specified");

	vector<string> errors(nThreads);

	if (nThreads == 1)
		startSingleThread(errors[0], numElements);
	else
	{
		vector<thread> threads(nThreads);
		vector<size_t> elementsInThread = getElementsPerThread(nThreads, numElements);

		for (size_t i = 0 ; i < threads.size() ; i++)
		{
			threads[i] = createOneThread(errors[i], elementsInThread[i]);
			updatePointers(elementsInThread[i]);
		}
		
		for (auto &t : threads)
			t.join();
	}

	for (auto &e : errors)
		if (e.length() > 0)
			return E(e);

	return true;
}

bool threadsTraceTheta(GravitationalLens &l, string &errStr,
		               double Ds, double Dds,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *betaX, double *betaY, size_t betaStride,
					   size_t numElements,
					   size_t nThreads)
{
	auto trace = [&l,Ds,Dds](string *errStr,
							   const double *thetaX, const double *thetaY, size_t thetaStride,
	                           double *betaX, double *betaY, size_t betaStride, size_t numElements)
	{
		for (size_t i = 0 ; i < numElements ; i++)
		{
			Vector2Dd t(*thetaX, *thetaY);
			Vector2Dd b;
			if (!l.traceTheta(Ds, Dds, t, &b))
			{
				*errStr = "Can't trace theta: " + l.getErrorString();
				break;
			}

			*betaX = b.getX();
			*betaY = b.getY();
			thetaX += thetaStride;
			thetaY += thetaStride;
			betaX += betaStride;
			betaY += betaStride;
		}
	};

	auto calcSingleThread = [&trace, thetaX, thetaY, thetaStride, betaX, betaY, betaStride](string &errStr, size_t numElements)
	{
		trace(&errStr, thetaX, thetaY, thetaStride, betaX, betaY, betaStride, numElements);
	};
	
	auto startNewThread = [&trace, &thetaX, &thetaY, thetaStride, &betaX, &betaY, betaStride](string &errStr, size_t numElements)
	{
		return thread(trace, &errStr, thetaX, thetaY, thetaStride, betaX, betaY, betaStride, numElements);
	};

	auto update = [&thetaX, &thetaY, thetaStride, &betaX, &betaY, betaStride](size_t numElements)
	{
		thetaX += numElements*thetaStride;
		thetaY += numElements*thetaStride;
		betaX += numElements*betaStride;
		betaY += numElements*betaStride;
	};

	return runThreads(l, errStr, numElements, nThreads, calcSingleThread, startNewThread, update);
}

bool threadsGetAlphaVector(grale::GravitationalLens &l, std::string &errStr,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *alphaX, double *alphaY, size_t alphaStride,
					   size_t numElements,
					   size_t nThreads)
{
	auto getAlpha = [&l](string *errStr,
							   const double *thetaX, const double *thetaY, size_t thetaStride,
	                           double *alphaX, double *alphaY, size_t alphaStride, size_t numElements)
	{
		for (size_t i = 0 ; i < numElements ; i++)
		{
			Vector2Dd t(*thetaX, *thetaY);
			Vector2Dd a;
			if (!l.getAlphaVector(t, &a))
			{
				*errStr = "Can't calculate alpha: " + l.getErrorString();
				break;
			}

			*alphaX = a.getX();
			*alphaY = a.getY();
			thetaX += thetaStride;
			thetaY += thetaStride;
			alphaX += alphaStride;
			alphaY += alphaStride;
		}
	};

	auto calcSingleThread = [&getAlpha, thetaX, thetaY, thetaStride, alphaX, alphaY, alphaStride](string &errStr, size_t numElements)
	{
		getAlpha(&errStr, thetaX, thetaY, thetaStride, alphaX, alphaY, alphaStride, numElements);
	};
	
	auto startNewThread = [&getAlpha, &thetaX, &thetaY, thetaStride, &alphaX, &alphaY, alphaStride](string &errStr, size_t numElements)
	{
		return thread(getAlpha, &errStr, thetaX, thetaY, thetaStride, alphaX, alphaY, alphaStride, numElements);
	};

	auto update = [&thetaX, &thetaY, thetaStride, &alphaX, &alphaY, alphaStride](size_t numElements)
	{
		thetaX += numElements*thetaStride;
		thetaY += numElements*thetaStride;
		alphaX += numElements*alphaStride;
		alphaY += numElements*alphaStride;
	};

	return runThreads(l, errStr, numElements, nThreads, calcSingleThread, startNewThread, update);
}
		
bool threadsGetAlphaVectorDerivatives(grale::GravitationalLens &l, std::string &errStr,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *alphaXX, double *alphaYY, double *alphaXY, size_t alphaStride,
					   size_t numElements,
					   size_t nThreads)
{
	auto getAlphaDerivs = [&l](string *errStr,
							   const double *thetaX, const double *thetaY, size_t thetaStride,
	                           double *alphaXX, double *alphaYY, double *alphaXY, 
							   size_t alphaStride, size_t numElements)
	{
		for (size_t i = 0 ; i < numElements ; i++)
		{
			Vector2Dd t(*thetaX, *thetaY);
			if (!l.getAlphaVectorDerivatives(t, *alphaXX, *alphaYY, *alphaXY))
			{
				*errStr = "Can't calculate alpha derivatives: " + l.getErrorString();
				break;
			}

			thetaX += thetaStride;
			thetaY += thetaStride;
			alphaXX += alphaStride;
			alphaYY += alphaStride;
			alphaXY += alphaStride;
		}
	};

	auto calcSingleThread = [&getAlphaDerivs, thetaX, thetaY, thetaStride, alphaXX, alphaYY, alphaXY, alphaStride](string &errStr, size_t numElements)
	{
		getAlphaDerivs(&errStr, thetaX, thetaY, thetaStride, alphaXX, alphaYY, alphaXY, alphaStride, numElements);
	};
	
	auto startNewThread = [&getAlphaDerivs, &thetaX, &thetaY, thetaStride, &alphaXX, &alphaYY, &alphaXY, alphaStride](string &errStr, size_t numElements)
	{
		return thread(getAlphaDerivs, &errStr, thetaX, thetaY, thetaStride, alphaXX, alphaYY, alphaXY, alphaStride, numElements);
	};

	auto update = [&thetaX, &thetaY, thetaStride, &alphaXX, &alphaYY, &alphaXY, alphaStride](size_t numElements)
	{
		thetaX += numElements*thetaStride;
		thetaY += numElements*thetaStride;
		alphaXX += numElements*alphaStride;
		alphaYY += numElements*alphaStride;
		alphaXY += numElements*alphaStride;
	};

	return runThreads(l, errStr, numElements, nThreads, calcSingleThread, startNewThread, update);
}

bool threadsGetInverseMagnification(grale::GravitationalLens &l, std::string &errStr,
		               double Ds, double Dds,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *invMag, size_t invMagStride,
					   size_t numElements,
					   size_t nThreads)
{
	auto getInv = [&l, Ds, Dds](string *errStr,
							   const double *thetaX, const double *thetaY, size_t thetaStride,
	                           double *invMag, size_t invMagStride, size_t numElements)
	{
		for (size_t i = 0 ; i < numElements ; i++)
		{
			Vector2Dd t(*thetaX, *thetaY);
			*invMag = l.getInverseMagnification(Ds, Dds, t);

			thetaX += thetaStride;
			thetaY += thetaStride;
			invMag += invMagStride;
		}
	};

	auto calcSingleThread = [&getInv, thetaX, thetaY, thetaStride, invMag, invMagStride](string &errStr, size_t numElements)
	{
		getInv(&errStr, thetaX, thetaY, thetaStride, invMag, invMagStride, numElements);
	};
	
	auto startNewThread = [&getInv, &thetaX, &thetaY, thetaStride, &invMag, invMagStride](string &errStr, size_t numElements)
	{
		return thread(getInv, &errStr, thetaX, thetaY, thetaStride, invMag, invMagStride, numElements);
	};

	auto update = [&thetaX, &thetaY, thetaStride, &invMag, invMagStride](size_t numElements)
	{
		thetaX += numElements*thetaStride;
		thetaY += numElements*thetaStride;
		invMag += numElements*invMagStride;
	};

	return runThreads(l, errStr, numElements, nThreads, calcSingleThread, startNewThread, update);
}

bool threadsGetSurfaceMassDensity(grale::GravitationalLens &l, std::string &errStr,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *dens, size_t densStride,
					   size_t numElements,
					   size_t nThreads)
{
	auto getDens = [&l](string *errStr,
							   const double *thetaX, const double *thetaY, size_t thetaStride,
	                           double *dens, size_t densStride, size_t numElements)
	{
		for (size_t i = 0 ; i < numElements ; i++)
		{
			Vector2Dd t(*thetaX, *thetaY);
			*dens = l.getSurfaceMassDensity(t);

			thetaX += thetaStride;
			thetaY += thetaStride;
			dens += densStride;
		}
	};

	auto calcSingleThread = [&getDens, thetaX, thetaY, thetaStride, dens, densStride](string &errStr, size_t numElements)
	{
		getDens(&errStr, thetaX, thetaY, thetaStride, dens, densStride, numElements);
	};
	
	auto startNewThread = [&getDens, &thetaX, &thetaY, thetaStride, &dens, densStride](string &errStr, size_t numElements)
	{
		return thread(getDens, &errStr, thetaX, thetaY, thetaStride, dens, densStride, numElements);
	};

	auto update = [&thetaX, &thetaY, thetaStride, &dens, densStride](size_t numElements)
	{
		thetaX += numElements*thetaStride;
		thetaY += numElements*thetaStride;
		dens += numElements*densStride;
	};

	return runThreads(l, errStr, numElements, nThreads, calcSingleThread, startNewThread, update);
}

bool threadsGetProjectedPotential(grale::GravitationalLens &l, std::string &errStr,
		               double Ds, double Dds,
					   const double *thetaX, const double *thetaY, size_t thetaStride,
					   double *potential, size_t potentialStride,
					   size_t numElements,
					   size_t nThreads)
{
	auto getPsi = [&l, Ds, Dds](string *errStr,
							   const double *thetaX, const double *thetaY, size_t thetaStride,
	                           double *potential, size_t potentialStride, size_t numElements)
	{
		for (size_t i = 0 ; i < numElements ; i++)
		{
			Vector2Dd t(*thetaX, *thetaY);
			if (!l.getProjectedPotential(Ds, Dds, t, potential))
			{
				*errStr = "Can't calculate potential: " + l.getErrorString();
				break;
			}

			thetaX += thetaStride;
			thetaY += thetaStride;
			potential += potentialStride;
		}
	};

	auto calcSingleThread = [&getPsi, thetaX, thetaY, thetaStride, potential, potentialStride](string &errStr, size_t numElements)
	{
		getPsi(&errStr, thetaX, thetaY, thetaStride, potential, potentialStride, numElements);
	};
	
	auto startNewThread = [&getPsi, &thetaX, &thetaY, thetaStride, &potential, potentialStride](string &errStr, size_t numElements)
	{
		return thread(getPsi, &errStr, thetaX, thetaY, thetaStride, potential, potentialStride, numElements);
	};

	auto update = [&thetaX, &thetaY, thetaStride, &potential, potentialStride](size_t numElements)
	{
		thetaX += numElements*thetaStride;
		thetaY += numElements*thetaStride;
		potential += numElements*potentialStride;
	};

	return runThreads(l, errStr, numElements, nThreads, calcSingleThread, startNewThread, update);
}

