#include <mpi.h>
#include "communicator.h"
#include "inputoutput.h"
#include "gravitationallens.h"
#include <serut/memoryserializer.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>

using namespace std;
using namespace grale;
using namespace serut;
using namespace errut;

class MPIRenderer : public Communicator
{
public:
	MPIRenderer() { }
	~MPIRenderer() { }

	static void runMPIHelper();
protected:
	string getRenderType() const { return "LENSPLANE"; }
	string getVersionInfo() const { return "MPI mass density renderer"; }
private:
	bool_t renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
	                  double x0, double y0, double dX, double dY, int numX, int numY,
	                  vector<double> &renderPoints);
	bool_t renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
	                         const std::vector<double> &inputXY, std::vector<double> &renderPoints);

	void setProgress(double current, double target)												{ Communicator::setProgress(current, target); }
	static void getPixelOffsets(int numXY, vector<int> &offsets);
	static void getOffsetsAndSendCounts(int numXY, vector<int> &offsets, vector<int> &sendCounts);
	static void renderNodeResults(MPIRenderer *pComm, GravitationalLens *pLens, double x0, double y0, 
			                        double dX, double dY,
                                    int numX, int numY, vector<double> &nodeResults, vector<int> &offsets,
									double *pTarget);
	static void renderNodeResults(MPIRenderer *pComm, GravitationalLens *pLens, int numXY, const double *pXYPoints,
									double *pTarget);
};

void MPIRenderer::getPixelOffsets(int numXY, vector<int> &offsets)
{
	int worldSize;

	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	offsets.resize(worldSize+1);

	int total = numXY;
	double pixelsPerNode = (double)total/(double)worldSize;
	double off = pixelsPerNode;

	offsets[0] = 0;
	for (int i = 1 ; i < worldSize ; i++, off += pixelsPerNode)
		offsets[i] = (int)(off + 0.5);

	offsets[worldSize] = total;
}

void MPIRenderer::renderNodeResults(MPIRenderer *pComm, GravitationalLens *pLens, double x0, double y0, double dX, double dY,
                                    int numX, int numY, vector<double> &nodeResults, vector<int> &offsets,
									double *pTarget)
{
	getPixelOffsets(numX*numY, offsets);
	//cerr << "Offsets: ";
	//for (int i = 0 ; i < offsets.size() ; i++)
	//	cerr << " " << offsets[i];
	//cerr << endl;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	assert(rank + 1 < offsets.size());
	int numPixels = offsets[rank+1] - offsets[rank];

	nodeResults.resize(numPixels);

	int reportInterval = (int)((double)numPixels/100.0 + 0.5);
	if (reportInterval < 10) // process at least 10 pixels before reporting
		reportInterval = 10;
	else if (reportInterval > numX) // report at most one line before reporting
		reportInterval = numX;

	for (int i = 0, j = offsets[rank] ; i < numPixels ; i++, j++)
	{
		int xi = j % numX;
		int yi = j / numX;

		double y = y0 + dY*yi;
		double x = x0 + dX*xi;

		Vector2Dd theta(x, y);

		nodeResults[i] = pLens->getSurfaceMassDensity(theta);
	
		if (pComm && i%reportInterval == 0)
			pComm->setProgress(i+1, numPixels);
	}

	if (pComm)
		pComm->setProgress(numPixels, numPixels);

	vector<int> recvCounts(offsets.size()-1);
	for (size_t i = 0 ; i < recvCounts.size() ; i++)
		recvCounts[i] = (offsets[i+1] - offsets[i]);

	//cerr << recvCounts[rank] << "," << nodeResults.size() << endl;

	MPI_Gatherv(&(nodeResults[0]), nodeResults.size(), MPI_DOUBLE, pTarget, 
			    &(recvCounts[0]), &(offsets[0]), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

bool_t MPIRenderer::renderGrid(const vector<uint8_t> &lensData, GravitationalLens *pLens, 
		                      double x0, double y0, double dX, double dY, int numX, int numY,
							  vector<double> &renderPoints)
{
	renderPoints.resize(numX*numY);

	// Inform helpers that this is the grid based method
	int renderingGrid = 1;
	MPI_Bcast(&renderingGrid, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Let's first distribute the lens
	int lensSize = (int)lensData.size();
	uint8_t *pLensData = (uint8_t*)&(lensData[0]);

	MPI_Bcast(&lensSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(pLensData, lensSize, MPI_BYTE, 0, MPI_COMM_WORLD);

	// Distribute parameters
	double dParams[4] = { x0, y0, dX, dY };
	int iParams[2] = { numX, numY };

	MPI_Bcast(dParams, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(iParams, 2, MPI_INT, 0, MPI_COMM_WORLD);

	// Let each node render it's parts
	vector<int> offsets;
	vector<double> nodeResults;
	renderNodeResults(this, pLens, x0, y0, dX, dY, numX, numY, nodeResults, offsets, &(renderPoints[0]));

	return true;
}

bool_t MPIRenderer::renderPointVector(const std::vector<uint8_t> &lensData, grale::GravitationalLens *pLens, 
		                              const std::vector<double> &inputXY, std::vector<double> &renderPoints)
{
	int numXY = inputXY.size()/2;
	renderPoints.resize(numXY);

	// Inform helpers that this is NOT the grid based method
	int renderingGrid = 0;
	MPI_Bcast(&renderingGrid, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Let's distribute the lens
	int lensSize = (int)lensData.size();
	uint8_t *pLensData = (uint8_t*)&(lensData[0]);

	MPI_Bcast(&lensSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(pLensData, lensSize, MPI_BYTE, 0, MPI_COMM_WORLD);

	MPI_Bcast(&numXY, 1, MPI_INT, 0, MPI_COMM_WORLD);

	renderNodeResults(this, pLens, numXY, &(inputXY[0]), &(renderPoints[0]));
	return true;
}

void MPIRenderer::getOffsetsAndSendCounts(int numXY, vector<int> &offsets, vector<int> &sendCounts)
{
	getPixelOffsets(numXY, offsets);
	for (auto &offset : offsets) // each pixel is two coords
		offset *= 2;

	sendCounts.resize(offsets.size()-1); // offsets is made one longer than number of processes
	for (int i = 0 ; i < offsets.size()-1 ; i++)
		sendCounts[i] = offsets[i+1] - offsets[i];
}

void MPIRenderer::renderNodeResults(MPIRenderer *pComm, GravitationalLens *pLens, int numXY, const double *pXYPoints,
									double *pTarget)
{
	vector<int> offsets, sendCounts;
	getOffsetsAndSendCounts(numXY, offsets, sendCounts); // note that offsets are multiplied by two to use as number of doubles

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numRecv = sendCounts[rank];
	vector<double> xyRenderPart(numRecv);

	MPI_Scatterv(pXYPoints, &sendCounts[0], &offsets[0], MPI_DOUBLE, &xyRenderPart[0], numRecv, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int pixelsToRender = numRecv/2;
	assert(numRecv%2 == 0);

	int reportInterval = (int)((double)pixelsToRender/100.0 + 0.5);
	if (reportInterval < 10) // process at least 10 pixels before reporting
		reportInterval = 10;

	vector<double> nodeResults(pixelsToRender);

	for (int i = 0 ; i < pixelsToRender ;  i++)
	{
		double x = xyRenderPart[i*2+0];
		double y = xyRenderPart[i*2+1];

		Vector2Dd theta(x, y);
		nodeResults[i] = pLens->getSurfaceMassDensity(theta);

		if (pComm && i%reportInterval == 0)
			pComm->setProgress(i+1, pixelsToRender);
	}

	if (pComm)
		pComm->setProgress(pixelsToRender, pixelsToRender);

	vector<int> doubleOffsets(offsets.size());

	// Offsets are just pixel offsets, but each node creates 5 doubles per pixel
	for (size_t i = 0 ; i < doubleOffsets.size() ; i++)
		doubleOffsets[i] = (offsets[i]/2); // dividing by two again because offsets had become number of doubles

	vector<int> recvCounts(offsets.size()-1);
	for (size_t i = 0 ; i < recvCounts.size() ; i++)
		recvCounts[i] = (doubleOffsets[i+1] - doubleOffsets[i]);

	//cerr << recvCounts[rank] << "," << nodeResults.size() << endl;

	MPI_Gatherv(&(nodeResults[0]), nodeResults.size(), MPI_DOUBLE, pTarget, 
			    &(recvCounts[0]), &(doubleOffsets[0]), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MPIRenderer::runMPIHelper()
{
	int renderingGrid = 0;
	MPI_Bcast(&renderingGrid, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Read the lens data
	int lensSize = 0;
	MPI_Bcast(&lensSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	vector<uint8_t> lensData(lensSize);
	uint8_t *pLensData = (uint8_t*)&(lensData[0]);

	MPI_Bcast(pLensData, lensSize, MPI_BYTE, 0, MPI_COMM_WORLD);

	GravitationalLens *pLens = 0;
	MemorySerializer mSer(&(lensData[0]), lensSize, 0, 0);
	string errStr;

	if (!GravitationalLens::read(mSer, &pLens, errStr))
	{
		cerr << errStr << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	//cerr << "Lens loaded successfully in helper" << endl;
	
	if (renderingGrid)
	{
		double dParams[4];
		int iParams[2];

		MPI_Bcast(dParams, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(iParams, 2, MPI_INT, 0, MPI_COMM_WORLD);

		// Let each node render it's parts
		vector<int> offsets;
		vector<double> nodeResults;
		renderNodeResults(0, pLens, dParams[0], dParams[1], dParams[2], dParams[3], iParams[0], iParams[1], nodeResults, offsets, NULL);
	}
	else
	{
		int numXY = 0;
		MPI_Bcast(&numXY, 1, MPI_INT, 0, MPI_COMM_WORLD);

		renderNodeResults(0, pLens, numXY, 0, 0);
	}
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	int worldSize, myRank;

	MPI_Comm_size( MPI_COMM_WORLD, &worldSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	if (argc != 3)
	{
		cerr << "ERROR: the names of input and output pipes are needed on command line" << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	if (myRank == 0) // The root node will communicate with the python script
	{
		string inputPipeName(argv[1]);
		string outputPipeName(argv[2]);
		cerr << "Opening " << outputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdOutFileDescriptor = open(outputPipeName.c_str(), O_WRONLY); 
		cerr << "Opening " << inputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = open(inputPipeName.c_str(), O_RDWR); 

		cerr << "MPI worldSize = " << worldSize << endl;

		MPIRenderer renderer;
		bool_t r;
		
		if (!(r = renderer.render()))
			cerr << "ERROR: " << r.getErrorString() << endl;
	}
	else
	{
		MPIRenderer::runMPIHelper();
	}

	MPI_Finalize();
	return 0;
}
