#pragma once

#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#endif // !WIN32
#include <memory>

template<class ParentClass>
class ThreadCommunicator : public ParentClass
{
public:
	ThreadCommunicator(size_t numThreads) : m_numThreads(numThreads)
	{
		std::cerr << "Using " << numThreads << " threads " << std::endl;
	}

	~ThreadCommunicator() { }
protected:
	std::string getVersionInfo() const override { return "EATk Thread based algorithm, " + std::to_string(m_numThreads) + " threads"; }

	errut::bool_t getCalculator(const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									eatk::IndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) override
	{
		return ParentClass::getMultiThreadedPopulationCalculator(m_numThreads, lensFitnessObjectType, calculatorType,
				                                    calcFactory, genomeCalculator, factoryParamBytes, calc);
	}
private:
	size_t m_numThreads;
};

template <class ParentClass>
int main_template(int argc, char *argv[])
{
	grale::LOG.init(argv[0]);
#ifndef WIN32
	// This was added to be able to launch instance in gdb, in which case
	// stdin/stdout can no longer be used for communication with the
	// executable (needed to control gdb).
	if (argc == 3)
	{
		std::string inputPipeName(argv[1]);
		std::string outputPipeName(argv[2]);
		std::cerr << "Opening " << outputPipeName << " for communication" << std::endl;
		// Nasty: override the file desc for communication
		stdOutFileDescriptor = open(outputPipeName.c_str(), O_WRONLY); 
		std::cerr << "Opening " << inputPipeName << " for communication" << std::endl;
		// Nasty: override the file desc for communication
		stdInFileDescriptor = open(inputPipeName.c_str(), O_RDWR); 
	}
#endif // WIN32

	size_t numThreads = 1;
	if (std::getenv("GRALE_NUMTHREADS"))
		numThreads = std::stoi(getenv("GRALE_NUMTHREADS"));

	ThreadCommunicator<ParentClass> comm(numThreads);

	errut::bool_t r = comm.run();
	if (!r)
	{
		std::cerr << "ERROR: " << r.getErrorString() << std::endl;
		return -1;
	}

	return 0;
}


