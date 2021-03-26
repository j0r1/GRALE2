#ifdef WIN32
#include <io.h>
#include <fcntl.h>
#endif
#include "log.h"
#include "inversioncommunicatornewga.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactorycommon.h"
#include "lensinversiongenome.h"
#include "galensmodule.h"
#include "gravitationallens.h"
#include "utils.h"
#include "lensgaindividual.h"
#include <serut/memoryserializer.h>
#include <serut/vectorserializer.h>
#include <serut/dummyserializer.h>
#include <mogal/geneticalgorithm.h>
#include <mogal/gamodule.h>
#include <memory>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace serut;
using namespace errut;
using namespace grale;

InversionCommunicator::InversionCommunicator()
{
#ifdef WIN32
	_setmode(_fileno(stdin), _O_BINARY);
	_setmode(_fileno(stdout), _O_BINARY);
#endif // WIN32
	m_nds = false;
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
	string moduleDir;
	if (!getenv("GRALE2_MODULEPATH", moduleDir))
		return "Environment variable GRALE2_MODULEPATH is not set, don't know where to look for inversion modules";

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

	LOG(Log::DBG, "Waiting for EXIT line");
	string exitStr;
	if (!(r = ReadLineStdin(10000, exitStr)))
		return "Couldn't read 'EXIT' line: " + r.getErrorString();

	LOG(Log::DBG, "Read '" + exitStr + "'");
	if (exitStr != "EXIT")
		return "Expected 'EXIT', but got: '" + exitStr + "'";

	LOG(Log::DBG, "At end of run");
	return true;
}

bool_t InversionCommunicator::runModule(const string &moduleDir, const string &moduleFile, GALensModule *pModule)
{
	bool_t r;
	int popSize = 0;

	auto pBaseFactory = pModule->createFactoryInstance();
	if (!pBaseFactory)
		return "Unable to create GA factory instance";
	unique_ptr<mogal::GAFactory> baseFactory(pBaseFactory); // just to make memory management easier

	auto pFactory = dynamic_cast<LensInversionGAFactoryCommon*>(pBaseFactory);
	if (!pFactory)
		return "GA factory instance could be created from module, but is not of expected type";

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
	mogal::GAFactoryParams *pFactoryParams = pFactory->createParamsInstance();
	if (!pFactoryParams)
		return "Unable to create factory parameters instance: " + pFactory->getErrorString();

	unique_ptr<mogal::GAFactoryParams> factoryParamsAutoDelete(pFactoryParams);

	if (!(r = loadFromBytes(gaParams, gaParamBytes)))
		return "Unable to load GA parameters from received data: " + r.getErrorString();

	if (!(r = loadFromBytes(*pFactoryParams, factoryParamBytes)))
		return "Unable to load GA factory parameters from received data: " + r.getErrorString();

	string runStr;
	if (!(r = ReadLineStdin(10000, runStr)))
		return "Unable to read 'RUN' command: " + r.getErrorString();
	if (runStr == "RUN")
		m_nds = false;
	else if (runStr == "RUN:NDS")
		m_nds = true;
	else
		return "Expected 'RUN' or 'RUN:NDS', but got: '" + runStr + "'";

	if (!pFactory->init(pFactoryParams))
		return "Unable to initialize the GA factory: " + pFactory->getErrorString();

	if (!(r = runGAWrapper(popSize, *pFactory, gaParams, moduleDir, moduleFile, *pModule, factoryParamBytes)))
		return r;

	LOG(Log::DBG, "Finished runGA");

	if (hasGeneticAlgorithm())
	{
		if (!(r = onGAFinished(*pFactory)))
			return r;
	}
	LOG(Log::DBG, "Called 'onGAFinished', end of runModule");
	return true;
}

bool_t InversionCommunicator::runGAWrapper(int popSize, mogal::GAFactory &factory, grale::GAParameters &params,
	                     const string &moduleDir, const string &moduleFile, grale::GALensModule &module,
						 const vector<uint8_t> &factoryParamBytes)
{
	bool_t r = runGA(popSize, factory, params, moduleDir, moduleFile, module, factoryParamBytes);
	if (!r)
	{
		string factoryError = factory.getErrorString();
		if (factoryError.length() > 0)
			return r.getErrorString() + "(factory has error: " + factoryError + ")";
		return r;
	}
	return true;
}

bool_t InversionCommunicator::runGA(int popSize, mogal::GAFactory &factory, grale::GAParameters &params,
	                     const string &moduleDir, const string &moduleFile, grale::GALensModule &module,
						 const vector<uint8_t> &factoryParamBytes)
{
	return "Not implemented in base class";
}

bool_t InversionCommunicator::onGAFinished(const LensInversionGAFactoryCommon &factory)
{
	vector<shared_ptr<mogal2::Individual>> bestGenomes;

	if (!m_nds)
	{
		auto g = getPreferredBestGenome();
		if (!g.get())
			return "No preferred best genome was available!";
		bestGenomes.push_back(g);
	}
	else // will send back entire non-dominated set
	{
		getAllBestGenomes(bestGenomes);
		if(bestGenomes.size() == 0)
			return "No non-dominated set wat stored!";
	}

	string fitnessCriteria = "TODO";
	
	WriteLineStdout("CRITERIA:" + fitnessCriteria);
	WriteLineStdout(strprintf("NUMSOLS:%d", (int)bestGenomes.size()));
	for (size_t i = 0 ; i < bestGenomes.size() ; i++)
	{
		grale::LensGAIndividual *pInd = dynamic_cast<grale::LensGAIndividual*>(bestGenomes[i].get());
		grale::LensGAGenome *pGenome = dynamic_cast<grale::LensGAGenome*>(pInd->genomePtr());
		if (!pInd || !pGenome)
			return "A genome in the best genomes set is not of expected type";

		string fitnessValues = pInd->fitness()->toString();
		LOG(Log::DBG, "Selected genome has fitness: " + fitnessValues);

		string errStr;

		LOG(Log::DBG, "Getting lens for genome");
		unique_ptr<GravitationalLens> lens(factory.createLens(pGenome->m_weights, pGenome->m_sheets, pGenome->m_scaleFactor, errStr));
		if (!lens.get())
			return "Unable to create lens from genome: " + errStr;

		LOG(Log::DBG, "Writing lens to buffer");
		VectorSerializer bufSer;
		if (!lens->write(bufSer))
			return "Error serializing resulting lens: " + lens->getErrorString();

		LOG(Log::DBG, "Writing buffer over connection");
		auto buf = bufSer.getBuffer();
		WriteLineStdout(strprintf("RESULT:%d", (int)buf.size()));
		WriteBytesStdout(buf);

		LOG(Log::DBG, "Writing fitness");
		WriteLineStdout("FITNESS:" + fitnessValues);
	}

	LOG(Log::DBG, "End of onGAFinished");
	return true;
}
