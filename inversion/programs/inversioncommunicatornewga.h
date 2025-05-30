#pragma once

#include "graleconfig.h"
#include "lensgagenomecalculator.h"
#include "lensfitnessobject.h"
#include "lensgacalculatorregistry.h"
#include "lensgamultipopulationparameters.h"
#include <errut/booltype.h>
#include <serut/memoryserializer.h>
#include <vector>
#include <stdint.h>
#include <eatk/population.h>

namespace grale
{
	class LensInversionGAFactoryCommon;
	class GALensModule;
	class EAParameters;
	class LensGAConvergenceParameters;
}

class InversionCommunicator
{
public:
	typedef errut::bool_t bool_t;

	InversionCommunicator();
	virtual ~InversionCommunicator();

	bool_t run();

	template<class T> static  bool_t loadFromBytes(T &x, const std::vector<uint8_t> &bytes);
protected:
	virtual std::string getVersionInfo() const = 0;
	virtual bool_t runModule(const std::string &lensFitnessObjectType, 
	                         std::unique_ptr<grale::LensFitnessObject> fitnessObject,
							 const std::string &calculatorType);

	virtual bool_t runGA(int popSize, const std::string &lensFitnessObjectType, 
						 const std::string &calculatorType,
	                     grale::LensGACalculatorFactory &calcFactory, 
						 const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
						 const std::vector<uint8_t> &factoryParamBytes,
						 const std::vector<std::unique_ptr<grale::EAParameters>> &allEAParams,
						 const std::vector<grale::LensGAConvergenceParameters> &convParams,
						 const std::shared_ptr<grale::LensGAMultiPopulationParameters> &multiPopParams,
						 const std::vector<std::string> &allEATypes);

	bool_t readLineWithPrefix(const std::string &prefix, std::string &value, int timeoutMSec);
	bool_t readLineWithPrefix(const std::string &prefix, int &value, int timeoutMSec);
	bool_t readLineAndBytesWithPrefix(const std::string &prefix, std::vector<uint8_t> &bytes, int timeoutMSec);

	virtual bool hasGeneticAlgorithm() const { return false; }
	virtual void getAllBestGenomes(std::vector<std::shared_ptr<eatk::Individual>> &bestGenomes) { }
	virtual std::shared_ptr<eatk::Individual> getPreferredBestGenome() { return nullptr; }

	bool_t onGAFinished(const grale::LensGAGenomeCalculator &calculator);
	bool m_nds;
};

template<class T>
inline InversionCommunicator::bool_t InversionCommunicator::loadFromBytes(T &x, const std::vector<uint8_t> &bytes)
{
	serut::MemorySerializer mSer(&bytes[0], bytes.size(), 0, 0);
	if (!x.read(mSer))
		return x.getErrorString();
	return true;
}
