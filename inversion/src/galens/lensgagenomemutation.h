#pragma once

#include "graleconfig.h"
#include <eatk/crossovermutation.h>
#include <eatk/randomnumbergenerator.h>

namespace grale
{

class LensGAGenomeMutation  : public eatk::GenomeMutation
{
public:
	LensGAGenomeMutation(const std::shared_ptr<eatk::RandomNumberGenerator> &rng, float chanceMultiplier,
					   bool allowNegativeValues, float mutationAmplitude, bool absoluteMutation);

	void setMutationAmplitude(float x) { m_mutationAmplitude = x; }
	void setAbsoluteMutation(bool x) { m_absoluteMutation = x; }

	errut::bool_t check(const eatk::Genome &genome) override;
	errut::bool_t mutate(eatk::Genome &genome, bool &isChanged) override;
private: 
	std::shared_ptr<eatk::RandomNumberGenerator> m_rng;
	float m_chanceMultiplier;
	float m_mutationAmplitude;
	bool m_allowNegative;
	bool m_absoluteMutation;
};

}