#pragma once

#include "graleconfig.h"
#include <mogal2/crossovermutation.h>
#include <mogal2/randomnumbergenerator.h>

namespace grale
{

class LensGAGenomeMutation  : public mogal2::GenomeMutation
{
public:
	LensGAGenomeMutation(const std::shared_ptr<mogal2::RandomNumberGenerator> &rng, float chanceMultiplier,
					   bool allowNegativeValues, float mutationAmplitude, bool absoluteMutation);

	void setMutationAmplitude(float x) { m_mutationAmplitude = x; }
	void setAbsoluteMutation(bool x) { m_absoluteMutation = x; }

	errut::bool_t check(const mogal2::Genome &genome) override;
	errut::bool_t mutate(mogal2::Genome &genome, bool &isChanged) override;
private: 
	std::shared_ptr<mogal2::RandomNumberGenerator> m_rng;
	float m_chanceMultiplier;
	float m_mutationAmplitude;
	bool m_allowNegative;
	bool m_absoluteMutation;
};

}