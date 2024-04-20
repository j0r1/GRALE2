#pragma once

#include <eatk/migrationstrategy.h>

class MyExchange : public eatk::SequentialRandomIndividualExchange
{
public:
	MyExchange(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, size_t iterations) : eatk::SequentialRandomIndividualExchange(rng, iterations) { }
protected:
	void onExchange(size_t generation, size_t srcPop, size_t srcIndividualIdx, size_t dstPop, size_t dstIndividualIdx) override
	{
		std::cerr << "Generation " << generation << ": migrating " << srcIndividualIdx << " from pop " << srcPop << " to " << dstIndividualIdx << " in pop " << dstPop << std::endl;
	}
};


