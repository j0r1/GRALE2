#ifdef WIN32
#include <io.h>
#include <fcntl.h>
#endif
#include "inversioncommunicator.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "gridlensinversiongafactoryparams.h"
#include "galensmodule.h"
#include "lensinversiongenomebase.h"
#include "gravitationallens.h"
#include <serut/memoryserializer.h>
#include <serut/dummyserializer.h>
#include <mogal/geneticalgorithm.h>
#include <mogal/gamodule.h>
#include <memory>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace grale;
using namespace serut;
using namespace mogal;
using namespace errut;

InversionCommunicator::InversionCommunicator()
{
#ifdef WIN32
	_setmode(_fileno(stdin), _O_BINARY);
	_setmode(_fileno(stdout), _O_BINARY);
#endif // WIN32
}

InversionCommunicator::~InversionCommunicator()
{
}

bool_t InversionCommunicator::readLineWithPrefix(const string &prefix, string &value, int timeoutMSec)
{
	string line;
	string prefix2 = prefix + ":";
	bool_t r;

	if (!(r = ReadLineStdin(timeoutMSec, line)))
		return "Error reading line with prefix '" + prefix + "': " + r.getErrorString();

	if (!startsWith(line, prefix2))
		return "The read line does not start with prefix " + prefix2;
	value = line.substr(prefix2.length());
	return true;
}

bool_t InversionCommunicator::readLineWithPrefix(const string &prefix, int &value, int timeoutMSec)
{
	string valueStr;
	bool_t r = readLineWithPrefix(prefix, valueStr, timeoutMSec);
	if (!r)
		return r;

	if (!parseAsInt(valueStr, value))
		return "Unable to interpret '" + valueStr + "' as an integer number";
	return true;
}

bool_t InversionCommunicator::readLineAndBytesWithPrefix(const string &prefix, vector<uint8_t> &bytes, int timeoutMSec)
{
	bool_t r;
	int numBytes = 0;

	if (!(r = readLineWithPrefix(prefix, numBytes, timeoutMSec)))
		return r;
	
	if (numBytes < 0)
		return "A negative number of bytes was specified in reading " + prefix;
	if (numBytes > 0)
	{
		bytes.resize(numBytes);
		if (!(r = ReadBytesStdin(bytes)))
			return r;
	}

	return true;
}

bool_t InversionCommunicator::run()
{
	const char envName[] = "GRALE2_MODULEPATH";
	if (getenv(envName) == nullptr)
		return "Environment variable GRALE2_MODULEPATH is not set, don't know where to look for inversion modules";

	string moduleDir(getenv(envName));
	bool_t r;

	if (!(r = WriteLineStdout("GAINVERTER:" + getVersionInfo())))
		return "Unable to send identification: " + r.getErrorString();

	string moduleFile;
	if (!(r = readLineWithPrefix("MODULE", moduleFile, 60000)))
		return "Error reading module name: " + r.getErrorString();

	GALensModule module;
	if (!module.open(moduleDir, moduleFile))
		return "Unable to open " + moduleFile + " in directory " + moduleDir + ": " + module.getErrorString();

	if (!(r = runModule(moduleDir, moduleFile, &module)))
		return r;

	string exitStr;
	if (!(r = ReadLineStdin(10000, exitStr)))
		return "Couldn't read 'EXIT' line: " + r.getErrorString();
	if (exitStr != "EXIT")
		return "Expected 'EXIT', but got: '" + exitStr + "'";

	return true;
}

bool_t InversionCommunicator::runModule(const string &moduleDir, const string &moduleFile, GALensModule *pModule)
{
	bool_t r;
	int popSize = 0;

	if (!(r = readLineWithPrefix("POPULATIONSIZE", popSize, 10000)))
		return "Error reading population size: " + r.getErrorString();
	if (popSize < 1)
		return "A negative population size was specified";

	vector<uint8_t> gaParamBytes;
	if (!(r = readLineAndBytesWithPrefix("GAPARAMS", gaParamBytes, 10000)))
		return "Error reading GA parameters: " + r.getErrorString();
	
	vector<uint8_t> factoryParamBytes;
	if (!(r = readLineAndBytesWithPrefix("GAFACTORYPARAMS", factoryParamBytes, 10000)))
		return "Error reading GA factory parameters: " + r.getErrorString();

	GAParameters gaParams;
	GridLensInversionGAFactoryParams factoryParams;
	if (!(r = loadFromBytes(gaParams, gaParamBytes)))
		return "Unable to load GA parameters from received data: " + r.getErrorString();
	if (!(r = loadFromBytes(factoryParams, factoryParamBytes)))
		return "Unable to load GA factory parameters from received data: " + r.getErrorString();

	string runStr;
	if (!(r = ReadLineStdin(10000, runStr)))
		return "Unable to read 'RUN' command: " + r.getErrorString();
	if (runStr != "RUN")
		return "Expected 'RUN', but got: '" + runStr + "'";

	GeneticAlgorithmParams params(gaParams.getSelectionPressure(), gaParams.getUseElitism(),
	                              gaParams.getAlwaysIncludeBest(), gaParams.getCrossOverRate());

	auto pFactory = pModule->createFactoryInstance();
	if (!pFactory)
		return "Unable to create GA factory instance";

	unique_ptr<GAFactory> factory(pFactory); // just to make memory management easier

	if (!pFactory->init(&factoryParams))
		return "Unable to initialize the GA factory: " + pFactory->getErrorString();

	if (!(r = runGA(popSize, *pFactory, params, moduleDir, moduleFile, factoryParamBytes)))
		return r;

	return true;
}

bool_t InversionCommunicator::runGA(int popSize, GAFactory &factory, GeneticAlgorithmParams &params,
	                     const string &moduleDir, const string &moduleFile,
						 const vector<uint8_t> &factoryParamBytes)
{
	return "Not implemented in base class";
}

bool_t InversionCommunicator::onGAFinished(mogal::GeneticAlgorithm &ga)
{
	const LensInversionGenomeBase *pGenome = static_cast<const LensInversionGenomeBase *>(ga.selectPreferredGenome());
	if (!pGenome)
		return "No best genome was set: " + ga.getErrorString();
	
	string fitnessCriteria = "TODO";
	string fitnessValues = pGenome->getFitnessDescription();
	cout << "Selected genome has fitness: "	<< fitnessValues << endl;
	cout.flush();
	
	string errStr;
	double totalmass; // This isn't really used anymore

	unique_ptr<GravitationalLens> lens(pGenome->createLens(&totalmass, errStr));
	if (!lens.get())
		return "Unable to create lens from genome: " + errStr;

	DummySerializer dSer;
	if (!lens->write(dSer))
		return "Error serializing resulting lens: " + lens->getErrorString();

	vector<uint8_t> buf(dSer.getBytesWritten());

	MemorySerializer mSer(nullptr, 0, &buf[0], buf.size());
	if (!lens->write(mSer))
		return "Error serializing resulting lens: " + lens->getErrorString();

	// Write the result to stdout
	WriteLineStdout(strprintf("RESULT:%d", (int)buf.size()));
	WriteBytesStdout(buf);
	WriteLineStdout("FITNESS:" + fitnessValues);
	WriteLineStdout("CRITERIA:" + fitnessCriteria);

	return true;
}

