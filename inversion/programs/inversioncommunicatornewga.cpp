#ifdef WIN32
#include <io.h>
#include <fcntl.h>
#endif
#include "log.h"
#include "inversioncommunicatornewga.h"
#include "inputoutput.h"
#include "gaparameters.h"
#include "lensinversiongafactorycommon.h"
#include "gravitationallens.h"
#include "utils.h"
#include "lensgaindividual.h"
#include "lensfitnessobject.h"
#include <serut/memoryserializer.h>
#include <serut/vectorserializer.h>
#include <serut/dummyserializer.h>
#include <memory>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace serut;
using namespace errut;
using namespace grale;

class MyCalcLogger : public grale::LensGAGenomeCalculatorLogger
{
public:
	void log(const std::string &s) const override
	{
		WriteLineStdout("GAMESSAGESTR:" + s);
	}
};


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
		return "The read line does not start with prefix " + prefix2 + " (" + line + ")";
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
	// string moduleDir;
	// if (!getenv("GRALE2_MODULEPATH", moduleDir))
	// 	return "Environment variable GRALE2_MODULEPATH is not set, don't know where to look for inversion modules";

	bool_t r;
	if (!(r = WriteLineStdout("GAINVERTER:" + getVersionInfo())))
		return "Unable to send identification: " + r.getErrorString();

	string lensFitnessObjectType;
	if (!(r = readLineWithPrefix("FITNESSOBJECT", lensFitnessObjectType, 60000)))
		return "Error reading fitnessObject name: " + r.getErrorString();

	string calculatorType;
	if (!(r = readLineWithPrefix("CALCULATOR", calculatorType, 60000)))
		return "Error reading calculator type: " + r.getErrorString();

	unique_ptr<LensFitnessObject> fitObj = LensFitnessObjectRegistry::instance().createFitnessObject(lensFitnessObjectType);
	if (!fitObj.get())
		return "No fitness object with name '" + lensFitnessObjectType + "' is known";

	if (!(r = runModule(lensFitnessObjectType, move(fitObj), calculatorType)))
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

bool_t InversionCommunicator::runModule(const std::string &lensFitnessObjectType,
										std::unique_ptr<grale::LensFitnessObject> fitnessObject,
                                        const std::string &calculatorType)
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
	
	vector<uint8_t> calcParamBytes;
	if (!(r = readLineAndBytesWithPrefix("GAFACTORYPARAMS", calcParamBytes, 10000)))
		return "Error reading GA factory parameters: " + r.getErrorString();

	string runStr;
	if (!(r = ReadLineStdin(10000, runStr)))
		return "Unable to read 'RUN' command: " + r.getErrorString();
	if (runStr == "RUN")
		m_nds = false;
	else if (runStr == "RUN:NDS")
		m_nds = true;
	else
		return "Expected 'RUN' or 'RUN:NDS', but got: '" + runStr + "'";

	LensGACalculatorFactory *pCalculatorFactory = LensGACalculatorRegistry::instance().getFactory(calculatorType);
	if (!pCalculatorFactory)
		return "Calculator type '" + calculatorType + "' is not known";

	GAParameters gaParams;
	if (!(r = loadFromBytes(gaParams, gaParamBytes)))
		return "Unable to load GA parameters from received data: " + r.getErrorString();

	auto calculatorParams = pCalculatorFactory->createParametersInstance();
	if (!(r = loadFromBytes(*calculatorParams, calcParamBytes)))
		return "Can't load calculator parameters from received data: " + r.getErrorString();

	shared_ptr<grale::LensGAGenomeCalculator> calculatorInstance = pCalculatorFactory->createCalculatorInstance(move(fitnessObject));
	auto logger = make_shared<MyCalcLogger>();
	calculatorInstance->setLogger(logger);

	if (!(r = calculatorInstance->init(*calculatorParams)))
		return "Unable to initialize calculator: " + r.getErrorString();

	if (!(r = runGA(popSize, lensFitnessObjectType, calculatorType, *pCalculatorFactory,
	                calculatorInstance, calcParamBytes, gaParams)))
		return r;

	LOG(Log::DBG, "Finished runGA");

	if (hasGeneticAlgorithm())
	{
		if (!(r = onGAFinished(*calculatorInstance)))
			return r;
	}
	LOG(Log::DBG, "Called 'onGAFinished', end of runModule");
	return true;
}

bool_t InversionCommunicator::runGA(int popSize, const std::string &lensFitnessObjectType, 
	                     const std::string &calcType, grale::LensGACalculatorFactory &calcFactory,
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
						 const std::vector<uint8_t> &factoryParamBytes, const grale::GAParameters &params)
{
	return "Not implemented in base class";
}

bool_t InversionCommunicator::onGAFinished(const grale::LensGAGenomeCalculator &calculator)
{
	vector<shared_ptr<eatk::Individual>> bestGenomes;

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
		
		bool_t r;
		unique_ptr<grale::GravitationalLens> lens;
		if (!(r = calculator.createLens(*pGenome, lens)))
			return "Unable to create lens from genome: " + r.getErrorString();

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
