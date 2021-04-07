#include "newgacommunicatorbase.h"
#include "inversionregistry.h"

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
	string getVersionInfo() const override { return "MOGAL2 Thread based algorithm, " + to_string(m_numThreads) + " threads"; }

	bool_t getCalculator(const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									grale::LensGAIndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) override
	{
		bool_t r;
		vector<shared_ptr<eatk::GenomeFitnessCalculation>> genomeFitnessCalculators = { genomeCalculator };

		for (size_t i = 1 ; i < m_numThreads ; i++)
		{
			unique_ptr<grale::LensFitnessObject> fitObj = grale::LensFitnessObjectRegistry::instance().createFitnessObject(lensFitnessObjectType);
			if (!fitObj.get())
				return "No fitness object with name '" + lensFitnessObjectType + "' is known";

			auto calculatorParams = calcFactory.createParametersInstance();
			if (!(r = InversionCommunicator::loadFromBytes(*calculatorParams, factoryParamBytes)))
				return "Can't load calculator parameters from received data: " + r.getErrorString();
		
	 		shared_ptr<grale::LensGAGenomeCalculator> calculatorInstance = calcFactory.createCalculatorInstance(move(fitObj));
			if (!(r = calculatorInstance->init(*calculatorParams)))
			 	return "Unable to initialize calculator: " + r.getErrorString();

			genomeFitnessCalculators.push_back(calculatorInstance);
		}
		
		if (m_numThreads <= 1)
			calc = make_shared<eatk::SingleThreadedPopulationFitnessCalculation>(genomeFitnessCalculators[0]);
		else
		{
			auto mpCalc = make_shared<eatk::MultiThreadedPopulationFitnessCalculation>();
			if (!(r = mpCalc->initThreadPool(genomeFitnessCalculators)))
				return "Unable to initialize threads: " + r.getErrorString();
			calc = mpCalc;
		}
		return true;
	}
private:
	size_t m_numThreads;
};

int main(int argc, char *argv[])
{
	grale::registerDefaultInversionComponents();

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
	if (getenv("NUMTHREADS"))
		numThreads = stoi(getenv("NUMTHREADS"));

	NewGACommunicatorThreads comm(numThreads);

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
