#pragma once

#include "randomnumbergenerator.h"
#include <mogal2/randomnumbergenerator.h>
#include <stdlib.h>
#include <iostream>

// Wrap the GSL based RNG for now
class GslRNGWrapper : public mogal2::RandomNumberGenerator
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
private:
	grale::RandomNumberGenerator m_rng;
};

