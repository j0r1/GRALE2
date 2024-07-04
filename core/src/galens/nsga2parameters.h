#pragma once

#include "graleconfig.h"
#include "eaparameters.h"

namespace grale
{

// Parameters using default crossover/mutation operators
class NSGA2Parameters : public EAParameters
{
public:
	NSGA2Parameters(double smallMutationSize = -1); // Negative means absolute mutation);
	~NSGA2Parameters();
	// Negative or zero means absolute mutation
	double getSmallMutationSize() const													{ return m_smallMutSize; }

	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	double m_smallMutSize;
};

}