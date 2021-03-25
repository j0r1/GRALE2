#include <mpi.h>
#include "newgacommunicatorbase.h"
#include <serut/memoryserializer.h>
#include <fcntl.h>
#include <unistd.h>

#define MPIABORT(errMsg) do { grale::LOG(grale::Log::ERR, (errMsg)); cerr << errMsg << std::endl; cerr.flush(); sleep(10); MPI_Abort(MPI_COMM_WORLD, -1); } while(0)

using namespace std;
using namespace serut;
using namespace errut;

class NewGACommunicatorMPI : public NewGACommunicatorBase
{
public:
	NewGACommunicatorMPI(size_t s) : m_size(s) { }
	~NewGACommunicatorMPI() { }
protected:
	string getVersionInfo() const override { return "MOGAL2 MPI based algorithm, " + to_string(m_size) + " processes"; }

	bool_t getCalculator(grale::LensInversionGAFactoryCommon &gaFactory,
						const std::string &moduleDir, const std::string &moduleFile, grale::GALensModule &module,
						const std::vector<uint8_t> &factoryParamBytes,
						shared_ptr<mogal2::PopulationFitnessCalculation> &calc) override
	{
		bool_t r;

		// TODO
		return true;
	}
private:
	size_t m_size;
};

bool_t runHelper()
{
	int paramsSize = 0;
	MPI_Bcast(&paramsSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	vector<uint8_t> paramBuffer(paramsSize);
	MPI_Bcast(&(paramBuffer[0]), paramsSize, MPI_BYTE, 0, MPI_COMM_WORLD);
	MemorySerializer mSer(&(paramBuffer[0]), paramsSize, 0, 0);

	string baseDir, moduleName;
	
	if (!mSer.readString(baseDir))
		return "Internal error: can't read module base directory";
	if (!mSer.readString(moduleName))
		return "Internal error: can't read module name";
	
	mogal::GAModule module;

	if (!module.open(baseDir, moduleName))
		return "Internal error: can't open specified GA module: " + baseDir + "/" + moduleName + ": " + module.getErrorString();

	mogal::GAFactory *pFactory = module.createFactoryInstance();
	if (!pFactory)
		return "Internal error: unable to create factory from module: " + module.getErrorString();
	unique_ptr<mogal::GAFactory> factory(pFactory);

	mogal::GAFactoryParams *pParams = pFactory->createParamsInstance();
	if (!pParams)
		return "Internal error: unable to create factory parameters instance: " + pFactory->getErrorString();
	unique_ptr<mogal::GAFactoryParams> factParams(pParams);

	if (!pParams->read(mSer))
		return "Internal error: unable to read factory parameters: " + pParams->getErrorString();

	if (!pFactory->init(pParams))
		return "Internal error: unable to init factory: " + pFactory->getErrorString();

	// TODO

	return true;
}

int main(int argc, char *argv[])
{
	int worldSize, myRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	grale::LOG.init(argv[0]);
	
	if (myRank == 0)
		grale::LOG(grale::Log::INF, "Root node");
	else
		grale::LOG(grale::Log::INF, "Helper node");

	if (argc != 3)
		MPIABORT("ERROR: the names of input and output pipes are needed on command line");

	if (worldSize <= 1)
		MPIABORT("ERROR: Need more than one process to be able to work");

	if (myRank == 0)
	{
		string inputPipeName(argv[1]);
		string outputPipeName(argv[2]);
		cerr << "Opening " << outputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdOutFileDescriptor = open(outputPipeName.c_str(), O_WRONLY); 
		cerr << "Opening " << inputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = open(inputPipeName.c_str(), O_RDWR); 

		cerr << "Starting MPI session with " << worldSize << " processes" << endl;

		NewGACommunicatorMPI comm(worldSize);

		bool_t r = comm.run();
		if (!r)
			MPIABORT("ERROR: " + r.getErrorString());
	}
	else
	{
		bool_t r = runHelper();
		if (!r)
			MPIABORT("ERROR: " + r.getErrorString());
	}

	MPI_Finalize();
	return 0;
}
