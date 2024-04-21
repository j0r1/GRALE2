#pragma once

#include <eatk/randomnumbergenerator.h>
#include <memory>
#include <iostream>

class RngWrapper : public eatk::RandomNumberGenerator
{
public:
	RngWrapper(const std::shared_ptr<eatk::RandomNumberGenerator> &rng) : m_rng(rng) { }

	double getRandomDouble() override
	{
		double d = m_rng->getRandomDouble();
		std::cerr << "RNGWRAPPER:getRandomDouble() = " << d << std::endl;
		//if (std::abs(d-0.449974) < 1e-6)
		//	std::cerr << "JORI: first RNG correspondence after first phase" << std::endl;
		//if (std::abs(d-0.393911) < 1e-6)
		//	std::cerr << "JORI: last RNG correspondence" << std::endl;
		return d;
	}

	double getRandomDouble(double min, double max) override
	{
		double d = m_rng->getRandomDouble(min, max);
		std::cerr << "RNGWRAPPER:getRandomDouble(" << min << "," << max << ") = " << d << std::endl;
		return d;
	}

	float getRandomFloat() override
	{
		float f = m_rng->getRandomFloat();
		std::cerr << "RNGWRAPPER:getRandomFloat() = " << f << std::endl;
		//if (std::abs(f-0.27235) < 1e-6f)
		//	throw std::runtime_error("RngWrapper::getRandomFloat trigger");
		return f;
	}

	float getRandomFloat(float min, float max) override
	{
		float f = m_rng->getRandomFloat(min, max);
		std::cerr << "RNGWRAPPER:getRandomFloat(" << min << "," << max << ") = " << f << std::endl;
		return f;
	}

	uint32_t getRandomUint32() override
	{
		uint32_t x = m_rng->getRandomUint32();
		std::cerr << "RNGWRAPPER:getRandomUint32() = " << x << std::endl;
		return x;
	}
private:
	std::shared_ptr<eatk::RandomNumberGenerator> m_rng;
};
