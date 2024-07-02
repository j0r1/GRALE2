#pragma once

#include "graleconfig.h"
#include "lensgaindividual.h"
#include <eatk/crossovermutation.h>
#include <eatk/randomnumbergenerator.h>
#include <eatk/vectorgenomedelikecrossover.h>

namespace grale
{

class LensGAGenomeCrossover : public eatk::GenomeCrossover
{
public:
	LensGAGenomeCrossover(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, bool allowNegative);

	errut::bool_t check(const std::vector<std::shared_ptr<eatk::Genome>> &parents) override;
	errut::bool_t generateOffspring(const std::vector<std::shared_ptr<eatk::Genome>> &parents,
	                                std::vector<std::shared_ptr<eatk::Genome>> &generatedOffspring) override;
private:
	std::shared_ptr<eatk::RandomNumberGenerator> m_rng;
	bool m_allowNegative;
};

typedef eatk::DELikeCrossOverTemplate<float, LensGAGenome> LensGAGenomeDELikeCrossover;

}
