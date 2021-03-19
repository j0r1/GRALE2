#pragma once

#include "graleconfig.h"
#include <mogal2/crossovermutation.h>
#include <mogal2/randomnumbergenerator.h>

namespace grale
{

class LensGAGenomeCrossover : public mogal2::GenomeCrossover
{
public:
	LensGAGenomeCrossover(const std::shared_ptr<mogal2::RandomNumberGenerator> &rng, bool allowNegative);

	errut::bool_t check(const std::vector<std::shared_ptr<mogal2::Genome>> &parents) override;
	errut::bool_t generateOffspring(const std::vector<std::shared_ptr<mogal2::Genome>> &parents,
	                                std::vector<std::shared_ptr<mogal2::Genome>> &generatedOffspring) override;
private:
	std::shared_ptr<mogal2::RandomNumberGenerator> m_rng;
	bool m_allowNegative;
};

}
