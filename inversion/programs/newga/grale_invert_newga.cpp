#include "newgacommunicatorbase.h"
#include <serut/vectorserializer.h>
#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#endif // !WIN32

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

	bool_t getCalculator(grale::LensInversionGAFactoryCommon &gaFactory,
						const std::string &moduleDir, const std::string &moduleFile, grale::GALensModule &module,
						const std::vector<uint8_t> &factoryParamBytes,
						shared_ptr<mogal2::PopulationFitnessCalculation> &calc) override
	{
		bool_t r;
		vector<shared_ptr<mogal2::GenomeFitnessCalculation>> genomeFitnessCalculators = { make_shared<LensFitnessCalculation>(gaFactory) };

		for (size_t i = 2 ; i < m_numThreads ; i++)
		{
			auto newFac = shared_ptr<grale::LensInversionGAFactoryCommon>(dynamic_cast<grale::LensInversionGAFactoryCommon*>(module.createFactoryInstance()));
			auto facParams = shared_ptr<mogal::GAFactoryParams>(newFac->createParamsInstance());

			serut::VectorSerializer ser(factoryParamBytes);
			if (!facParams->read(ser))
				return "Coudln't deserialize factory parameters: " + facParams->getErrorString();

			if (!newFac->init(facParams.get()))
				return "Unable to init factory copy: " + newFac->getErrorString();
			
			genomeFitnessCalculators.push_back(make_shared<LensFitnessCalculation>(newFac));
		}
		
		if (m_numThreads <= 1)
			calc = make_shared<mogal2::SingleThreadedPopulationFitnessCalculation>(genomeFitnessCalculators[0]);
		else
		{
			auto mpCalc = make_shared<mogal2::MultiThreadedPopulationFitnessCalculation>();
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
