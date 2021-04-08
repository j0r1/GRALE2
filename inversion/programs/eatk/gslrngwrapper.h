#pragma once

#include "randomnumbergenerator.h"
#include <eatk/randomnumbergenerator.h>
#include <stdlib.h>
#include <iostream>

// Wrap the GSL based RNG for now
class GslRNGWrapper : public eatk::RandomNumberGenerator
{
public:
	GslRNGWrapper() { }
	~GslRNGWrapper() { }
    double getRandomDouble() override { return m_rng.pickRandomNumber(); }
    float getRandomFloat() override { return (float)m_rng.pickRandomNumber(); }
	uint32_t getRandomUint32() override
	{
		std::cerr << "getRandomUint32 NOT IMPLEMENTED" << std::endl;
		exit(-1);
		return 0;
	}
	uint32_t getSeed() const { return m_rng.getSeed(); }
private:
	grale::RandomNumberGenerator m_rng;
};

