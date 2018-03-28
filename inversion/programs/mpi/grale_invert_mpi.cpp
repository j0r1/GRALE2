#include <mogal/mogalconfig.h>

#ifdef MOGALCONFIG_MPISUPPORT

#include <mpi.h>
#include "inputoutput.h"
#include "gaparameters.h"
#include "gridlensinversiongafactoryparams.h"
#include "inversioncommunicator.h"
#include <errut/booltype.h>
#include <mogal/mpigeneticalgorithm.h>
#include <mogal/gamodule.h>
#include <serut/memoryserializer.h>
#include <serut/dummyserializer.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <memory>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;
using namespace grale;
using namespace serut;
using namespace mogal;
using namespace errut;

#include "commonga.h"
typedef CommonGA<MPIGeneticAlgorithm> MyGA;

class MultiCoreCommunicator : public InversionCommunicator
{
public:
	MultiCoreCommunicator(int worldSize) { m_worldSize = worldSize; }
	~MultiCoreCommunicator() { }
protected:
	string getVersionInfo() const override 
	{
		return strprintf("Multi core (MPI) genetic algorithm engine (%d processes)", m_worldSize);
	}

	bool_t runGA(int popSize, GAFactory &factory, GeneticAlgorithmParams &params,
	             const std::string &moduleDir, const std::string &moduleFile,
	             const std::vector<uint8_t> &factoryParamBytes) override
	{
		DummySerializer dumSer;
		dumSer.writeString(moduleDir);
		dumSer.writeString(moduleFile);
		dumSer.writeBytes(&factoryParamBytes[0], factoryParamBytes.size());

		int paramSize = dumSer.getBytesWritten();
		
		vector<uint8_t> buffer(paramSize);
		MemorySerializer mSer(0, 0, &(buffer[0]), buffer.size());
		mSer.writeString(moduleDir);
		mSer.writeString(moduleFile);
		mSer.writeBytes(&factoryParamBytes[0], factoryParamBytes.size());

		MPI_Bcast(&paramSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&(buffer[0]), paramSize, MPI_BYTE, 0, MPI_COMM_WORLD);

		MyGA ga;
		if (!ga.runMPI(factory, popSize, &params))
			return "Unable to run MPI based GA: " + ga.getErrorString();

		return onGAFinished(ga);
	}
private:
	int m_worldSize;
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
	
	GAModule module;

	if (!module.open(baseDir, moduleName))
		return "Internal error: can't open specified GA module: " + baseDir + "/" + moduleName + ": " + module.getErrorString();

	GAFactory *pFactory = module.createFactoryInstance();
	if (!pFactory)
		return "Internal error: unable to create factory from module: " + module.getErrorString();
	unique_ptr<GAFactory> factory(pFactory);

	GAFactoryParams *pParams = pFactory->createParamsInstance();
	if (!pParams)
		return "Internal error: unable to create factory parameters instance: " + pFactory->getErrorString();
	unique_ptr<GAFactoryParams> factParams(pParams);

	if (!pParams->read(mSer))
		return "Internal error: unable to read factory parameters: " + pParams->getErrorString();

	if (!pFactory->init(pParams))
		return "Internal error: unable to init factory: " + pFactory->getErrorString();

	MPIGeneticAlgorithm gaHelper;

	if (!gaHelper.runMPIHelper(*pFactory))
		return "Internal error: error running MPI GA helper: " + gaHelper.getErrorString();

	return true;
}

int main(int argc, char *argv[])
{
	int worldSize, myRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	if (argc != 3)
	{
		cerr << "ERROR: the names of input and output pipes are needed on command line" << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	if (worldSize <= 1)
	{
		cerr << "ERROR: Need more than one process to be able to work" << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

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

		MultiCoreCommunicator comm(worldSize);

		bool_t r = comm.run();
		if (!r)
		{
			cerr << "ERROR: " << r.getErrorString() << endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}
	else
	{
		bool_t r = runHelper();
		if (!r)
		{
			cerr << "ERROR: " + r.getErrorString() << endl;
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}

	return 0;
}

#else

#include <iostream>

using namespace std;

int main(void)
{
	cerr << "ERROR: Can't run MPI genetic algorithm: MPI support was not enabled when compiling MOGAL" << endl;
	return -1;
}

#endif // MOGALCONFIG_MPISUPPORT
