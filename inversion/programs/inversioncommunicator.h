#ifndef INVERSIONCOMMUNICATOR_H

#define INVERSIONCOMMUNICATOR_H

#include "graleconfig.h"
#include <errut/booltype.h>
#include <serut/memoryserializer.h>
#include <vector>
#include <stdint.h>

namespace grale
{
	class GridLensInversionGAFactoryParams;
	class GALensModule;
}
namespace mogal
{
	class GeneticAlgorithmParams;
	class GeneticAlgorithm;
	class GAFactory;
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
	virtual bool_t runModule(const std::string &moduleDir, const std::string &moduleFile, grale::GALensModule *pModule);
	virtual bool_t runGA(int popSize, mogal::GAFactory &factory, mogal::GeneticAlgorithmParams &params,
	                     const std::string &moduleDir, const std::string &moduleFile,
						 const std::vector<uint8_t> &factoryParamBytes);

	bool_t readLineWithPrefix(const std::string &prefix, std::string &value, int timeoutMSec);
	bool_t readLineWithPrefix(const std::string &prefix, int &value, int timeoutMSec);
	bool_t readLineAndBytesWithPrefix(const std::string &prefix, std::vector<uint8_t> &bytes, int timeoutMSec);
	template<class T> bool_t loadFromBytes(T &x, std::vector<uint8_t> &bytes);

	bool_t onGAFinished(mogal::GeneticAlgorithm &ga);
};

template<class T>
inline InversionCommunicator::bool_t InversionCommunicator::loadFromBytes(T &x, std::vector<uint8_t> &bytes)
{
	serut::MemorySerializer mSer(&bytes[0], bytes.size(), 0, 0);
	if (!x.read(mSer))
		return x.getErrorString();
	return true;
}

#endif // INVERSIONCOMMUNICATOR_H
