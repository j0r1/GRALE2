#include "newgacommunicatorbase.h"

#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#endif // !WIN32
#include <memory>

using namespace std;
using namespace serut;
using namespace errut;

class NewGACommunicatorThreads : public NewGACommunicatorBase
{
public:
	NewGACommunicatorThreads(size_t numThreads) : m_numThreads(numThreads)
	{
		cerr << "Using " << numThreads << " threads " << endl;
	}

	~NewGACommunicatorThreads() { }
protected:
	string getVersionInfo() const override { return "EATk Thread based algorithm, " + to_string(m_numThreads) + " threads"; }

	bool_t getCalculator(const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									grale::LensGAIndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) override
	{
		return getMultiThreadedPopulationCalculator(m_numThreads, lensFitnessObjectType, calculatorType,
				                                    calcFactory, genomeCalculator, factoryParamBytes, calc);
	}
private:
	size_t m_numThreads;
};

int main(int argc, char *argv[])
{
	grale::LOG.init(argv[0]);
#ifndef WIN32
	// This was added to be able to launch instance in gdb, in which case
	// stdin/stdout can no longer be used for communication with the
	// executable (needed to control gdb).
	if (argc == 3)
	{
		string inputPipeName(argv[1]);
		string outputPipeName(argv[2]);
		cerr << "Opening " << outputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdOutFileDescriptor = open(outputPipeName.c_str(), O_WRONLY); 
		cerr << "Opening " << inputPipeName << " for communication" << endl;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = open(inputPipeName.c_str(), O_RDWR); 
	}
#endif // WIN32

	size_t numThreads = 1;
	if (getenv("GRALE_NUMTHREADS"))
		numThreads = stoi(getenv("GRALE_NUMTHREADS"));

	NewGACommunicatorThreads comm(numThreads);

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
