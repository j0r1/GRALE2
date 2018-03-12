#include "inputoutput.h"
#include "gaparameters.h"
#include "gridlensinversiongafactoryparams.h"
#include "inversioncommunicator.h"
#include <errut/booltype.h>
#include <serut/memoryserializer.h>
#include <mogal/geneticalgorithm.h>

#include <iostream>
#include <string>

using namespace std;
using namespace grale;
using namespace serut;
using namespace mogal;
using namespace errut;

#include "commonga.h"
typedef CommonGA<GeneticAlgorithm> MyGA;

class SingleCoreCommunicator : public InversionCommunicator
{
public:
	SingleCoreCommunicator() { }
	~SingleCoreCommunicator() { }
protected:
	string getVersionInfo() const override { return "Single core genetic algorithm engine"; }

	bool_t runGA(int popSize, GAFactory &factory, GeneticAlgorithmParams &params,
	             const std::string &moduleDir, const std::string &moduleFile,
	             const std::vector<uint8_t> &factoryParamBytes) override
	{
		MyGA ga;

		if (!ga.run(factory, popSize,&params))
			return "Error running GA: " + ga.getErrorString();

		return onGAFinished(ga);
	}
};

int main(int argc, char *argv[])
{
	SingleCoreCommunicator comm;

	bool_t r = comm.run();
	if (!r)
	{
		cerr << "ERROR: " << r.getErrorString() << endl;
		return -1;
	}

	return 0;
}
