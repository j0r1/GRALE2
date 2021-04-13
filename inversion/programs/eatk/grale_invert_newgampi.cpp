#include <mpi.h>
#include "newgacommunicatorbase.h"
#include <eatk/mpieventdistributor.h>
#include <eatk/mpipopulationfitnesscalculation.h>
#include <eatk/singlethreadedpopulationfitnesscalculation.h>
#include <serut/memoryserializer.h>
#include <serut/vectorserializer.h>
#include <sys/socket.h>
#include <sys/un.h>
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
	string getVersionInfo() const override { return "EATk MPI based algorithm, " + to_string(m_size) + " processes"; }

	void calculatorCleanup() override
	{
		m_evtDist->signal(eatk::MPIEventHandler::Done);
		m_evtDist = nullptr;
	}

	bool_t getCalculator(const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									grale::LensGAIndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) override
	{
		bool_t r;
		VectorSerializer ser;

		if (!(ser.writeString(lensFitnessObjectType) && ser.writeString(calculatorType) &&
		      ser.writeBytes(factoryParamBytes.data(), factoryParamBytes.size())))
			return "Error writing module info and parameters: " + ser.getErrorString();

		int paramSize = ser.getBufferSize();
		MPI_Bcast(&paramSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast((void*)ser.getBufferPointer(), ser.getBufferSize(), MPI_BYTE, 0, MPI_COMM_WORLD);
		
		auto localCalc = make_shared<eatk::SingleThreadedPopulationFitnessCalculation>(genomeCalculator);
		
		// The event distribution mechanism uses a weak pointer, make sure this doesn't
		// get deallocated too soon by storing it in the class instance
		m_evtDist = make_shared<eatk::MPIEventDistributor>();
		auto mpiCalc = make_shared<eatk::MPIPopulationFitnessCalculation>(m_evtDist);
		calc = mpiCalc;

		auto refGenome = creation.createUnInitializedGenome();
		auto refFitness = creation.createEmptyFitness();

		if (!(r = mpiCalc->init(*refGenome, *refFitness, localCalc)))
		{
			m_evtDist = nullptr;
			return "Error initializing MPI calculator: " + r.getErrorString();
		}

		return true;
	}
private:
	size_t m_size;
	std::shared_ptr<eatk::MPIEventDistributor> m_evtDist;
};

bool_t runHelper()
{
	int paramsSize = 0;
	MPI_Bcast(&paramsSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	vector<uint8_t> paramBuffer(paramsSize);
	MPI_Bcast(&(paramBuffer[0]), paramsSize, MPI_BYTE, 0, MPI_COMM_WORLD);
	MemorySerializer mSer(&(paramBuffer[0]), paramsSize, 0, 0);

	string lensFitnessObjectType, calculatorType;
	
	if (!mSer.readString(lensFitnessObjectType))
		return "Internal error: can't read fitness object type";
	if (!mSer.readString(calculatorType))
		return "Internal error: can't read calculator type";
	
	unique_ptr<grale::LensFitnessObject> fitObj = grale::LensFitnessObjectRegistry::instance().createFitnessObject(lensFitnessObjectType);
	if (!fitObj.get())
		return "No fitness object with name '" + lensFitnessObjectType + "' is known";

	grale::LensGACalculatorFactory *pCalculatorFactory = grale::LensGACalculatorRegistry::instance().getFactory(calculatorType);
	if (!pCalculatorFactory)
		return "Calculator type '" + calculatorType + "' is not known";

	auto calculatorParams = pCalculatorFactory->createParametersInstance();
	if (!calculatorParams->read(mSer))
		return "Internal error: unable to read factory parameters: " + calculatorParams->getErrorString();

	bool_t r;
	shared_ptr<grale::LensGAGenomeCalculator> calculatorInstance = pCalculatorFactory->createCalculatorInstance(move(fitObj));
	if (!(r = calculatorInstance->init(*calculatorParams)))
		return "Unable to initialize calculator: " + r.getErrorString();

	auto localCalc = make_shared<eatk::SingleThreadedPopulationFitnessCalculation>(calculatorInstance);
	auto dist = make_shared<eatk::MPIEventDistributor>();
	auto calc = make_shared<eatk::MPIPopulationFitnessCalculation>(dist);

	dist->setHandler(eatk::MPIEventHandler::Calculation, calc);

	grale::LensGAGenome refGenome(0, 0); // will receive layout in 'init'
	grale::LensGAFitness refFitness(0); // will receive the number of components in 'init'

	if (!(r = calc->init(refGenome, refFitness, localCalc)))
		return "Error initializing MPI calculation: " + r.getErrorString();
	return dist->eventLoop();
}

bool_t getConnectedSocket(const string &sockName, int &s)
{
	cerr << "Opening " << sockName << " for communication" << endl;

	s = socket(AF_UNIX, SOCK_STREAM, 0);
	if (s < 0)
		return "Error creating UNIX socket: code " + to_string(s);

	struct sockaddr_un addr;
	memset(&addr, 0, sizeof(struct sockaddr_un));
	addr.sun_family = AF_UNIX;
	strncpy(addr.sun_path, sockName.c_str(), sizeof(addr.sun_path) - 1);
	addr.sun_path[sizeof(addr.sun_path) - 1] = 0;

	int status = connect(s, (struct sockaddr *)&addr, sizeof(addr));
	if (status < 0)
		return "Error connecting to " + string(addr.sun_path) + ": code " + to_string(status);

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

	if (argc != 2)
		MPIABORT("ERROR: the names of UNIX socket path is needed on command line");

	if (worldSize <= 1)
		MPIABORT("ERROR: Need more than one process to be able to work");

	if (myRank == 0)
	{
		int s;
		bool_t r = getConnectedSocket(argv[1], s);
		if (!r)
			MPIABORT("Can't get connected UNIX socket: " + r.getErrorString());

		// Nasty: override the file desc for communication
		stdOutFileDescriptor = s;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = s;

		cerr << "Starting MPI session with " << worldSize << " processes" << endl;

		NewGACommunicatorMPI comm(worldSize);

		if (!(r = comm.run()))
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
