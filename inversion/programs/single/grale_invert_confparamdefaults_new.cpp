#include "inversioncommunicatornewga.h"
#include "inputoutput.h"
#include "lensfitnessobject.h"
#include "configurationparameters.h"
#include <serut/memoryserializer.h>
#include <serut/dummyserializer.h>
#include <memory>
#include <iostream>
#include <string>

using namespace std;
using namespace grale;
using namespace serut;
using namespace errut;

class ConfDefaultsCommunicator : public InversionCommunicator
{
public:
	ConfDefaultsCommunicator() { }
	~ConfDefaultsCommunicator() { }
protected:
	string getVersionInfo() const override { return "Default configuration parameters for module"; }

	bool_t runModule(const std::string &lensFitnessObjectType, 
	                 std::unique_ptr<grale::LensFitnessObject> fitnessObject,
					 const std::string &calculatorType) override
    {
		auto &f = fitnessObject;
		if (!f.get())
			return "Unable to create lens fitness object for module: " + lensFitnessObjectType;

		unique_ptr<ConfigurationParameters> params(f->getDefaultParametersInstance());
		vector<uint8_t> paramsBytes;

		if (params.get())
		{
			DummySerializer dSer;
			if (!params->write(dSer))
				return "Error serializing the default parameters for the specified module: " + params->getErrorString();

			paramsBytes.resize(dSer.getBytesWritten());
			
			MemorySerializer mSer(nullptr, 0, &paramsBytes[0], paramsBytes.size());
			if (!params->write(mSer))
				return "Error serializing the default parameters for the specified module: " + params->getErrorString();
		}

		WriteLineStdout(strprintf("CONFPARAMDEFAULTS:%d", paramsBytes.size()));
		if (paramsBytes.size() > 0)
			WriteBytesStdout(&paramsBytes[0], paramsBytes.size());

		return true;
	}
};

int main(int argc, char *argv[])
{
	ConfDefaultsCommunicator comm;

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
