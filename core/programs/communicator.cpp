#define __STDC_FORMAT_MACROS
#ifdef WIN32
#include <io.h>
#include <fcntl.h>
#endif
#include "communicator.h"
#include "inputoutput.h"
#include "gravitationallens.h"
#include <serut/memoryserializer.h>
#include <inttypes.h>
#include <memory>
#include <iostream>

using namespace std;
using namespace serut;
using namespace grale;
using namespace errut;

Communicator::Communicator()
{
#ifdef WIN32
	_setmode(_fileno(stdin), _O_BINARY);
	_setmode(_fileno(stdout), _O_BINARY);
#endif // WIN32
}

Communicator::~Communicator()
{
}

bool_t Communicator::render()
{
	string line;
	bool_t r;

	if (!(r = WriteLineStdout("RENDERER:" + getVersionInfo() )))
		return "Unable to send identification: " + r.getErrorString();

	if (!(r = ReadLineStdin(10000, line)))
		return "Error reading requested renderer type: " + r.getErrorString();

	if (line != "LENSPLANE")
		return "Invalid renderer '" + line + "' requested";

	if (!(r = ReadLineStdin(1000, line)))
		return "Unable to read lens data: " + r.getErrorString();

	int lensSize = 0;
	vector<string> parts;
	SplitLine(line, parts);

	if (parts.size() != 2 || parts[0] != "LENSSIZE" || !parseAsInt(parts[1], lensSize))
		return "Expecting LENSSIZE numbytes, but got '" + line + "'";

	if (lensSize <= 0)
		return strprintf("Lens size should be positive, but was %d", lensSize);

	cerr << "Reading lens of " << lensSize << " bytes" << endl;
	vector<uint8_t> lensBytes(lensSize);
	if (!(r = ReadBytesStdin(lensBytes)))
		return "Error reading lens data: " + r.getErrorString();

	cerr << "Reading rendering type (grid/separate input points)" << endl;
	if (!(r = ReadLineStdin(1000, line)))
		return "Unable to read render type: " + r.getErrorString();

	int pointVectorSize = 0;

	SplitLine(line, parts);
	if (parts.size() != 3 || parts[0] != "RENDERTYPE" || !parseAsInt(parts[2], pointVectorSize))
		return "Expecting RENDERTYPE, but got '" + line + "'";

	bool useGridMethod = true;
	if (parts[1] == "GRID")
	{
		useGridMethod = true;
		if (pointVectorSize != 0)
			return "For GRID render type, the next parameter should be 0, but is '" + parts[2] + "'";
	}
	else if (parts[1] == "POINTVECTOR")
	{
		useGridMethod = false;
		if (pointVectorSize < 1)
			return "Illegal value of point vector size: '" + parts[2] + "'";
	}
	else
		return "Invalid rendertype, should be GRID or POINTVECTOR, but is '" + parts[1] + "'";

	GravitationalLens *pLens = 0;
	{
		MemorySerializer mSer(&(lensBytes[0]), lensSize, 0, 0);
		string errStr;

		if (!GravitationalLens::read(mSer, &pLens, errStr))
			return "Couldn't create lens from received data";
	}
	// Just to make sure the lens is automatically deleted
	unique_ptr<GravitationalLens> lensAutoDelete(pLens);
	vector<double> renderPoints;

	if (useGridMethod)
	{
		cerr << "Reading corner info and resolution" << endl;
		vector<uint8_t> renderInfo(sizeof(double)*4 + sizeof(int32_t)*2); // corners and resolution
		if (!(r = ReadBytesStdin(renderInfo)))
			return "Error reading corners and resolution: " + r.getErrorString();

		Vector2Dd bottomLeft, topRight;
		int32_t numX = 0, numY = 0;

		{
			MemorySerializer mSer(&(renderInfo[0]), renderInfo.size(), 0, 0);
			
			mSer.readDoubles(bottomLeft.getComponents(), 2);
			mSer.readDoubles(topRight.getComponents(), 2);
			mSer.readInt32(&numX);
			mSer.readInt32(&numY);
		}
		cerr << strprintf("Rendering %dx%d from (%g,%g) - (%g,%g)", numX, numY, bottomLeft.getX(), bottomLeft.getY(), topRight.getX(), topRight.getY()) << endl;
		if (numX <= 1 || numY <= 1)
			return strprintf("Invalid number of pixels (%d, %d) (must be at least two)", numX, numY);

		double dX = (topRight.getX()-bottomLeft.getX())/((double)(numX-1));
		double dY = (topRight.getY()-bottomLeft.getY())/((double)(numY-1));

		if (dX < 0 || dY < 0)
			return "Either the width or the height of the specified region is negative";

		setStatus("Rendering grid data");

		if (!(r = renderGrid(lensBytes, pLens, bottomLeft.getX(), bottomLeft.getY(), dX, dY, numX, numY, renderPoints)))
			return "Unable to render lens properties on specified grid: " + r.getErrorString();
	}
	else
	{
		vector<double> inputPoints(pointVectorSize*2); // *2 for each coordinate
		uint8_t *pBytes = (uint8_t*)(&inputPoints[0]);

		if (!(r = ReadBytesStdin(pBytes, inputPoints.size()*sizeof(double))))
			return "Error reading input points: " + r.getErrorString();

		setStatus("Rendering based on points vector");

		if (!(r = renderPointVector(lensBytes, pLens, inputPoints, renderPoints)))
			return "Unable to render lens properties on specified grid: " + r.getErrorString();
	}

	setStatus("Rendering finished");

	cerr << "Writing result" << endl;

	uint64_t numResultBytes = renderPoints.size()*sizeof(double);
	if (!(r = WriteLineStdout(strprintf("RESULT:%" PRIu64, numResultBytes))) ||
		!(r = WriteBytesStdout(&(renderPoints[0]), numResultBytes)) )
		return "Error communicating result: " + r.getErrorString();
	
	cerr << "Finishing, waiting for exit command" << endl;
	if (!(r = ReadLineStdin(60000, line)))
		return "Error waiting for EXIT command: " + r.getErrorString();

	if (line != "EXIT")
		return "Expecting 'EXIT' but got '" + line + "'";

	cerr << "Renderer exiting as expected." << endl;
	return true;
}

void Communicator::setStatus(const std::string &s)
{
	WriteLineStdout("STATUS:" + s);
}

void Communicator::setProgress(double current, double target)
{
	WriteLineStdout(strprintf("PROGRESS: %g %g", current, target));
}

