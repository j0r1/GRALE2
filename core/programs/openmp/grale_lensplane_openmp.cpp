#include <omp.h>
#include "communicator.h"
#include "gravitationallens.h"
#include "inputoutput.h"
#include <iostream>

using namespace std;
using namespace grale;
using namespace errut;

class OpenMPRenderer : public Communicator
{
public:
	OpenMPRenderer() { }
	~OpenMPRenderer() { }
protected:
	string getRenderType() const { return "LENSPLANE"; }
	string getVersionInfo() const { return "OpenMP lensplane renderer"; }
private:
	bool_t renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
	                  double x0, double y0, double dX, double dY, int numX, int numY,
	                  vector<double> &renderPoints);
	bool_t renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                          const std::vector<double> &inputXY, std::vector<double> &renderPoints);
};

bool_t OpenMPRenderer::renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
                                  double x0, double y0, double dX, double dY, int numX, int numY,
                                  vector<double> &renderPoints)
{
	cerr << "Number of threads: " << omp_get_max_threads() << endl;

	renderPoints.resize(numX*numY*5); // 5 components per point

	int rootCount = 0;

	#pragma omp parallel for schedule(static)
	for (int yi = 0 ; yi < numY ; yi++)
	{
		double y = y0 + dY*yi;

		for (int xi = 0 ; xi < numX ; xi++)
		{
			double x = x0 + dX*xi;
			Vector2Dd theta(x, y);
			Vector2Dd alpha;
			double axx = 0, ayy = 0, axy = 0;

			// TODO: error handling!
			//       can't just return from OpenMP like this
			/*
			if (!pLens->getAlphaVector(theta, &alpha) ||
				!pLens->getAlphaVectorDerivatives(theta, axx, ayy, axy))
				return strprintf("Can't calculate deflection or derivatives for point (%d, %d)", xi, yi);
			*/
			pLens->getAlphaVector(theta, &alpha);
			pLens->getAlphaVectorDerivatives(theta, axx, ayy, axy);
			
			int offset = (xi+yi*numX)*5;
			renderPoints[offset+0] = alpha.getX();
			renderPoints[offset+1] = alpha.getY();
			renderPoints[offset+2] = axx;
			renderPoints[offset+3] = ayy;
			renderPoints[offset+4] = axy;
		}

		if (omp_get_thread_num() == 0) // root thread only
		{
			rootCount++;
			setProgress(MIN(rootCount * omp_get_num_threads(), numY), numY);
		}
	}
	setProgress(numY, numY);

	return true;
}

bool_t OpenMPRenderer::renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                          const std::vector<double> &inputXY, std::vector<double> &renderPoints)
{
	cerr << "Number of threads: " << omp_get_max_threads() << endl;

	int numXY = inputXY.size()/2;
	renderPoints.resize(numXY*5); // 5 components per point

	int rootCount = 0;

	int pixelsToRenderApprox = numXY/omp_get_max_threads() + 1;
	int reportInterval = (int)((double)pixelsToRenderApprox/100.0 + 0.5);
	if (reportInterval < 10) // process at least 10 pixels before reporting
		reportInterval = 10;

	#pragma omp parallel for schedule(static)
	for (int i = 0 ; i < numXY ; i++)
	{
		double x = inputXY[i*2+0];
		double y = inputXY[i*2+1];

		Vector2Dd theta(x, y);
		Vector2Dd alpha;
		double axx = 0, ayy = 0, axy = 0;

		// TODO: error handling!
		//       can't just return from OpenMP like this
		/*
		if (!pLens->getAlphaVector(theta, &alpha) ||
			!pLens->getAlphaVectorDerivatives(theta, axx, ayy, axy))
			return strprintf("Can't calculate deflection or derivatives for point (%d, %d)", xi, yi);
		*/
		pLens->getAlphaVector(theta, &alpha);
		pLens->getAlphaVectorDerivatives(theta, axx, ayy, axy);
		
		int offset = i*5;
		renderPoints[offset+0] = alpha.getX();
		renderPoints[offset+1] = alpha.getY();
		renderPoints[offset+2] = axx;
		renderPoints[offset+3] = ayy;
		renderPoints[offset+4] = axy;

		if (omp_get_thread_num() == 0) // root thread only
		{
			rootCount++;
			if (rootCount%reportInterval == 0)
				setProgress(MIN(rootCount+1, pixelsToRenderApprox), pixelsToRenderApprox);
		}
	}
	setProgress(pixelsToRenderApprox, pixelsToRenderApprox);

	return true;
}

int main(int argc, char *argv[])
{
	OpenMPRenderer renderer;
	bool_t r;

	if (!(r = renderer.render()))
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}
	return 0;
}
