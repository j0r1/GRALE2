#pragma once

template<class ParentClass>
class ThreadCommunicator : public ParentClass
{
public:
	ThreadCommunicator(size_t numThreads) : m_numThreads(numThreads)
	{
		std::cerr << "Using " << numThreads << " threads " << std::endl;
	}

	~ThreadCommunicator() { }
protected:
	std::string getVersionInfo() const override { return "EATk Thread based algorithm, " + std::to_string(m_numThreads) + " threads"; }

	errut::bool_t getCalculator(const std::string &lensFitnessObjectType, const std::string &calculatorType,
									grale::LensGACalculatorFactory &calcFactory, 
									const std::shared_ptr<grale::LensGAGenomeCalculator> &genomeCalculator,
									const std::vector<uint8_t> &factoryParamBytes,
									eatk::IndividualCreation &creation,
									std::shared_ptr<eatk::PopulationFitnessCalculation> &calc) override
	{
		return ParentClass::getMultiThreadedPopulationCalculator(m_numThreads, lensFitnessObjectType, calculatorType,
				                                    calcFactory, genomeCalculator, factoryParamBytes, calc);
	}
private:
	size_t m_numThreads;
};

