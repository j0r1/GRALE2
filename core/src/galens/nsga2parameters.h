#pragma once

#include "graleconfig.h"
#include "eaparameters.h"
#include <limits>

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

// Parameters using experimental DE-like crossover operators
class NSGA2DELikeCrossoverParameters : public EAParameters
{
public:
	NSGA2DELikeCrossoverParameters(bool extraParent = true,
		float F = std::numeric_limits<float>::quiet_NaN(),
		float CR = std::numeric_limits<float>::quiet_NaN());
	~NSGA2DELikeCrossoverParameters();

	bool useExtraParent() const { return m_extraParent; }
	// NaN means that random number between 0 and 1 will be used
	float getF() const { return m_F; }
	float getCR() const { return m_CR; }

	errut::bool_t readInternal(serut::SerializationInterface &si) override;
	errut::bool_t writeInternal(serut::SerializationInterface &si) const override;
private:
	bool m_extraParent;
	float m_F, m_CR;
};

}