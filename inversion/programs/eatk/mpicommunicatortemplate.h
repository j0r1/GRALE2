#pragma once

#include <mpi.h>
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

template <class ParentClass>
class MPICommunicatorMPI : public ParentClass
{
public:
	MPICommunicatorMPI(size_t s, size_t numThreads) : m_size(s), m_numThreads(numThreads) { }
	~MPICommunicatorMPI() { }
protected:
	std::string getVersionInfo() const override { return "EATk MPI based algorithm, " + std::to_string(m_size) + " processes and " + std::to_string(m_numThreads) + " threads per process"; }

	void calculatorCleanup() override
	{
		if (m_evtDist)
			m_evtDist->signal(eatk::MPIEventHandler::Done);
		m_evtDist = nullptr;
	}

	errut::bool_t getCalculator(const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									eatk::IndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) override
	{
		errut::bool_t r;
		serut::VectorSerializer ser;

		if (!(ser.writeString(lensFitnessObjectType) && ser.writeString(calculatorType) &&
		      ser.writeBytes(factoryParamBytes.data(), factoryParamBytes.size())))
			return "Error writing module info and parameters: " + ser.getErrorString();

		int paramSize = ser.getBufferSize();
		MPI_Bcast(&paramSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast((void*)ser.getBufferPointer(), ser.getBufferSize(), MPI_BYTE, 0, MPI_COMM_WORLD);

		std::shared_ptr<eatk::PopulationFitnessCalculation> localCalc;
		if (!(r = ParentClass::getMultiThreadedPopulationCalculator(m_numThreads, lensFitnessObjectType, calculatorType,
						                               calcFactory, genomeCalculator, factoryParamBytes, localCalc)))
			return "Error creating calculator for MPI process: " + r.getErrorString();

		// The event distribution mechanism uses a weak pointer, make sure this doesn't
		// get deallocated too soon by storing it in the class instance
		m_evtDist = std::make_shared<eatk::MPIEventDistributor>();
		auto mpiCalc = std::make_shared<eatk::MPIPopulationFitnessCalculation>(m_evtDist);
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
	size_t m_size, m_numThreads;
	std::shared_ptr<eatk::MPIEventDistributor> m_evtDist;
};

template<class ParentClass>
inline errut::bool_t runHelper(size_t numThreads)
{
	int paramsSize = 0;
	MPI_Bcast(&paramsSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	vector<uint8_t> paramBuffer(paramsSize);
	MPI_Bcast(&(paramBuffer[0]), paramsSize, MPI_BYTE, 0, MPI_COMM_WORLD);
	MemorySerializer mSer(&(paramBuffer[0]), paramsSize, 0, 0);

	std::string lensFitnessObjectType, calculatorType;
	
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

	VectorSerializer calcParamsBuffer;
	if (!calculatorParams->write(calcParamsBuffer))
		return "Error writing calculator parameters to buffer in helper process" + calculatorParams->getErrorString();
	
	errut::bool_t r;
	shared_ptr<grale::LensGAGenomeCalculator> calculatorInstance = pCalculatorFactory->createCalculatorInstance(move(fitObj));
	if (!(r = calculatorInstance->init(*calculatorParams)))
		return "Unable to initialize calculator: " + r.getErrorString();

	std::shared_ptr<eatk::PopulationFitnessCalculation> localCalc;
	if (!(r = ParentClass::getMultiThreadedPopulationCalculator(numThreads, lensFitnessObjectType, calculatorType,
					                                                      *pCalculatorFactory, calculatorInstance, calcParamsBuffer.getBuffer(),
																		  localCalc)))
		return "Error creating calculator for MPI helper process: " + r.getErrorString();

	auto dist = std::make_shared<eatk::MPIEventDistributor>();
	auto calc = std::make_shared<eatk::MPIPopulationFitnessCalculation>(dist);

	dist->setHandler(eatk::MPIEventHandler::Calculation, calc);

	shared_ptr<eatk::Genome> refGenome = ParentClass::getReferenceGenomeForMPI(); // will receive layout in 'init'
	shared_ptr<eatk::Fitness> refFitness = ParentClass::getReferenceFitnessForMPI(); // will receive the number of components in 'init'

	if (!(r = calc->init(*refGenome, *refFitness, localCalc)))
		return "Error initializing MPI calculation: " + r.getErrorString();

	return dist->eventLoop();
}

inline errut::bool_t getConnectedSocket(const std::string &sockName, int &s)
{
	cerr << "Opening " << sockName << " for communication" << endl;

	s = socket(AF_UNIX, SOCK_STREAM, 0);
	if (s < 0)
		return "Error creating UNIX socket: code " + std::to_string(s);

	struct sockaddr_un addr;
	memset(&addr, 0, sizeof(struct sockaddr_un));
	addr.sun_family = AF_UNIX;
	strncpy(addr.sun_path, sockName.c_str(), sizeof(addr.sun_path) - 1);
	addr.sun_path[sizeof(addr.sun_path) - 1] = 0;

	int status = connect(s, (struct sockaddr *)&addr, sizeof(addr));
	if (status < 0)
		return "Error connecting to " + std::string(addr.sun_path) + ": code " + std::to_string(status);

	return true;
}

template<class ParentClass>
int main_template(int argc, char *argv[])
{
	int requested = MPI_THREAD_FUNNELED;
	int provided = -1;
	MPI_Init_thread(&argc, &argv, requested, &provided);
	if (provided != requested)
		MPIABORT("ERROR: Unable to get the requested threading capabilities for MPI");

	int worldSize, myRank;
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	grale::LOG.init(argv[0]);
	
	if (myRank == 0)
		grale::LOG(grale::Log::INF, "Root node");
	else
		grale::LOG(grale::Log::INF, "Helper node");

	if (argc != 3)
		MPIABORT("ERROR: the name of UNIX socket path is needed on command line as well as the number of threads");

	if (worldSize <= 1)
		MPIABORT("ERROR: Need more than one process to be able to work");

	int numThreads = stoi(argv[2]);
	if (numThreads <= 0)
		MPIABORT("ERROR: invalid number of thread (" + std::to_string(numThreads) + ") per MPI process");

	if (myRank == 0)
	{
		int s;
		errut::bool_t r = getConnectedSocket(argv[1], s);
		if (!r)
			MPIABORT("Can't get connected UNIX socket: " + r.getErrorString());

		// Nasty: override the file desc for communication
		stdOutFileDescriptor = s;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = s;

		cerr << "Starting MPI session with " << worldSize << " processes" << endl;

		MPICommunicatorMPI<ParentClass> comm(worldSize, (size_t)numThreads);

		if (!(r = comm.run()))
			MPIABORT("ERROR: " + r.getErrorString());
	}
	else
	{
		errut::bool_t r = runHelper<ParentClass>(numThreads);
		if (!r)
			MPIABORT("ERROR: " + r.getErrorString());
	}

	MPI_Finalize();
	return 0;
}
