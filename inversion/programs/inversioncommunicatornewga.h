#pragma once

#include "graleconfig.h"
#include "lensgagenomecalculator.h"
#include "lensfitnessobject.h"
#include "lensgacalculatorregistry.h"
#include <errut/booltype.h>
#include <serut/memoryserializer.h>
#include <vector>
#include <stdint.h>
#include <eatk/population.h>

namespace grale
{
	class LensInversionGAFactoryCommon;
	class GALensModule;
	class GAParameters;
	class LensGAConvergenceParameters;
}

class InversionCommunicator
{
public:
	typedef errut::bool_t bool_t;

	InversionCommunicator();
	virtual ~InversionCommunicator();

	bool_t run();
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
						 const grale::GAParameters &params,
						 const grale::LensGAConvergenceParameters &convParams);

	bool_t readLineWithPrefix(const std::string &prefix, std::string &value, int timeoutMSec);
	bool_t readLineWithPrefix(const std::string &prefix, int &value, int timeoutMSec);
	bool_t readLineAndBytesWithPrefix(const std::string &prefix, std::vector<uint8_t> &bytes, int timeoutMSec);
	template<class T> bool_t loadFromBytes(T &x, const std::vector<uint8_t> &bytes);

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
