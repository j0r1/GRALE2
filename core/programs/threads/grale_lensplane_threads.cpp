#include "communicator.h"
#include "gravitationallens.h"
#include "inputoutput.h"
#include "utils.h"
#include <iostream>
#include <memory>
#include <thread>

using namespace std;
using namespace grale;
using namespace errut;

class ThreadsRenderer : public Communicator
{
public:
	ThreadsRenderer(int numThreads) : m_numThreads(numThreads) { }
	~ThreadsRenderer() { }
protected:
	string getRenderType() const { return "LENSPLANE"; }
	string getVersionInfo() const { return "Threads lensplane renderer"; }
private:
	bool_t renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
	                  double x0, double y0, double dX, double dY, int numX, int numY,
	                  vector<double> &renderPoints);
	bool_t renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                          const std::vector<double> &inputXY, std::vector<double> &renderPoints);

	const int m_numThreads;
};

bool_t ThreadsRenderer::renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
                                  double x0, double y0, double dX, double dY, int numX, int numY,
                                  vector<double> &renderPoints)
{
	cerr << "Number of threads: " << m_numThreads << endl;

	renderPoints.resize(numX*numY*5); // 5 components per point

	// Make copies of the lens, to be very sure that we're being thread safe
	// (especially needed if GSL is used internally)
	const int nt = m_numThreads;
	vector<string> errMsg(nt);
	vector<bool> errors(nt, false);
	vector<GravitationalLens *> lenses;
	vector<shared_ptr<GravitationalLens>> allocatedLenses;

	lenses.push_back(pLens);
	for (int i = 1 ; i < nt ; i++)
	{
		unique_ptr<GravitationalLens> pLensCopy = pLens->createCopy();
		if (!pLensCopy.get())
			return "Unable to create copy of lens: " + pLens->getErrorString();
		
		allocatedLenses.push_back(move(pLensCopy));
		lenses.push_back(allocatedLenses.back().get());
	}

	int rootCount = 0;

	auto threadFunction = [&](int threadIdx)
	{
		for (int yi = threadIdx ; yi < numY ; yi += m_numThreads)
		{
			const int t = threadIdx;
			if (errors[t])
				continue;

			double y = y0 + dY*yi;

			for (int xi = 0 ; xi < numX ; xi++)
			{
				double x = x0 + dX*xi;
				Vector2Dd theta(x, y);
				Vector2Dd alpha;
				double axx = 0, ayy = 0, axy = 0;

				if (!lenses[t]->getAlphaVector(theta, &alpha) ||
					!lenses[t]->getAlphaVectorDerivatives(theta, axx, ayy, axy))
				{
					errMsg[t] = strprintf("Can't calculate deflection or derivatives for point (%d, %d)", xi, yi);
					errors[t] = true;
					break;
				}
				
				int offset = (xi+yi*numX)*5;
				renderPoints[offset+0] = alpha.getX();
				renderPoints[offset+1] = alpha.getY();
				renderPoints[offset+2] = axx;
				renderPoints[offset+3] = ayy;
				renderPoints[offset+4] = axy;
			}

			if (threadIdx == 0) // root thread only
			{
				rootCount++;
				setProgress(MIN(rootCount * m_numThreads, numY), numY);
			}
		}
	}; // end threadFunction lambda

	vector<thread> threads;
	for (int i = 0 ; i < m_numThreads ; i++)
		threads.push_back(thread(threadFunction, i));

	for (int i = 0 ; i < m_numThreads ; i++)
		threads[i].join();

	setProgress(numY, numY);

	for (int i = 1 ; i < nt ; i++)
	{
		if (errors[i])
			return errMsg[i];
	}

	return true;
}

bool_t ThreadsRenderer::renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                          const std::vector<double> &inputXY, std::vector<double> &renderPoints)
{
	cerr << "Number of threads: " << m_numThreads << endl;

	int numXY = inputXY.size()/2;
	renderPoints.resize(numXY*5); // 5 components per point

	int rootCount = 0;

	int pixelsToRenderApprox = numXY/m_numThreads + 1;
	int reportInterval = (int)((double)pixelsToRenderApprox/100.0 + 0.5);
	if (reportInterval < 10) // process at least 10 pixels before reporting
		reportInterval = 10;

	// Make copies of the lens, to be very sure that we're being thread safe
	// (especially needed if GSL is used internally)
	const int nt = m_numThreads;
	vector<string> errMsg(nt);
	vector<bool> errors(nt, false);
	vector<GravitationalLens *> lenses;
	vector<shared_ptr<GravitationalLens>> allocatedLenses;

	lenses.push_back(pLens);
	for (int i = 1 ; i < nt ; i++)
	{
		unique_ptr<GravitationalLens> pLensCopy = pLens->createCopy();
		if (!pLensCopy)
			return "Unable to create copy of lens: " + pLens->getErrorString();
		
		allocatedLenses.push_back(move(pLensCopy));
		lenses.push_back(allocatedLenses.back().get());
	}

	auto threadFunction = [&](int threadIdx)
	{
		for (int i = threadIdx ; i < numXY ; i += m_numThreads)
		{
			const int t = threadIdx;
			if (errors[t])
				continue;

			double x = inputXY[i*2+0];
			double y = inputXY[i*2+1];

			Vector2Dd theta(x, y);
			Vector2Dd alpha;
			double axx = 0, ayy = 0, axy = 0;

			if (!lenses[t]->getAlphaVector(theta, &alpha) ||
				!lenses[t]->getAlphaVectorDerivatives(theta, axx, ayy, axy))
			{
				errMsg[t] = strprintf("Can't calculate deflection or derivatives for point (%d)", i);
				errors[t] = true;
				continue;
			}
			
			int offset = i*5;
			renderPoints[offset+0] = alpha.getX();
			renderPoints[offset+1] = alpha.getY();
			renderPoints[offset+2] = axx;
			renderPoints[offset+3] = ayy;
			renderPoints[offset+4] = axy;

			if (threadIdx == 0) // root thread only
			{
				rootCount++;
				if (rootCount%reportInterval == 0)
					setProgress(MIN(rootCount+1, pixelsToRenderApprox), pixelsToRenderApprox);
			}
		}
	}; // end threadFunction lambda

	vector<thread> threads;
	for (int i = 0 ; i < m_numThreads ; i++)
		threads.push_back(thread(threadFunction, i));

	for (int i = 0 ; i < m_numThreads ; i++)
		threads[i].join();

	setProgress(pixelsToRenderApprox, pixelsToRenderApprox);

	for (int i = 1 ; i < nt ; i++)
	{
		if (errors[i])
			return errMsg[i];
	}

	return true;
}

int main(int argc, char *argv[])
{
	int numThreads;
	bool_t r = getenv("GRALE_NUMTHREADS", numThreads, 1);

	if (!r)
	{
		cerr << "ERROR: Error getting number of threads from environment variable GRALE_NUMTHREADS: " << r.getErrorString() << endl;
		return -1;
	}

	ThreadsRenderer renderer(numThreads);

	if (!(r = renderer.render()))
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		cout << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}
	return 0;
}
